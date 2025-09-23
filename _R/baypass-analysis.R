

library(tidyverse)
library(poolfstat)
library(seqinr)

if (file.exists(".Rprofile")) source(".Rprofile")

# source("_R/baypass-utils.R")
# # These functions not necessary:
# rm(compute_genetic_offset, compute.local.scores, fmd.dist, geno2YN, plot.omega, simulate.baypass, simulate.PCcorrelated.covariate)


#'
#' # also move everything into a single folder:
#' for ((r=1; r<=3; r++)); do
#'     cd run${r}
#'     for ((i=1; i<=5; i++)); do
#'         cd baypass-beta_1.0-run${r}-sub${i}
#'         mv *.log *.out ../
#'         cd ..
#'         rm -r baypass-beta_1.0-run${r}-sub${i}
#'     done
#'     cd ..
#' done
#' ```
#'



# bp_runs <- map(1:3, \(r) {
#     concatenate_res(dir = sprintf("~/_data/_baypass/beta/run%i", r),
#                     anaprefix = sprintf("tany-beta_1.0-run%i-sub", r),
#                     anasep = "",
#                     extension = "",
#                     nsubsets = 5,
#                     snpdet_prefix = "~/_data/_baypass/beta/tany.snpdet.sub",
#                     retrieve_pi_xtx = FALSE,
#                     retrieve_bfis = TRUE,
#                     retrieve_c2 = FALSE) |>
#         set_names(c("contig", "pos", "bf_db")) |>
#         as_tibble() |>
#         mutate(bf = 10^(bf_db / 10),
#                run = factor(r, levels = 1:3)) |>
#         select(run, everything())
# }) |>
#     list_rbind()


# =============================================================================*
# Reading baypass output ----
# =============================================================================*

baypass_dir <- "/Volumes/Lucas512/__to-transfer/baypass"


snp_df <- paste0(baypass_dir, "/baypass-inputs/tany.snpdet") |>
    read_table(col_types = "cdcc", col_names = c("contig", "pos", "ref", "alt")) |>
    # Have to do it this way bc some are in scientific notation in file:
    mutate(pos = as.integer(pos)) |>
    select(contig, pos)


beta_df <- crossing(b = 1, r = 1:10) |>
    pmap(\(b, r) {
        fn <- sprintf(paste0("%s/beta/baypass-beta%i-run%i/tany-beta%i-",
                             "run%i_summary_betai.out.gz"),
                      baypass_dir, b, r, b, r)
        fn |>
            read_table(col_types = "iidddd", skip = 1,
                       col_names = c("cov", "mrk", "m_beta", "sd_beta",
                                     "pip", "bf_db")) |>
            mutate(run = r)

    }) |>
    list_rbind() |>
    select(run, everything())

bf_df <- beta_df |>
    group_by(mrk) |>
    summarize(bf_db = median(bf_db),
              bf = 10^(bf_db / 10)) |>
    (\(x) {
        # CHecks to make sure it's okay to bind with `snp_df`:
        stopifnot(identical(x$mrk, 1:nrow(x)))
        stopifnot(identical(nrow(x), nrow(snp_df)))
        return(x)
    })() |>
    bind_cols(snp_df) |>
    select(contig, pos, everything()) |>
    # To make sure they stay in order of the assembly:
    mutate(contig = factor(contig, levels = paste0("CONTIG", 1:45)))

# SNPs exceeding 30 deciban threshold
# (Gautier et al. 2018, doi: 10.1016/j.cub.2018.08.023)
candi_df <- bf_df |> filter(bf_db >= 30)







# =============================================================================*
# Matching loci to genes ----
# =============================================================================*


#' This removes transcript id that's an extension separated by
#' either the last '-' or '.'
#' This converts the transcript id to a gene id.
trans2genes <- function(trans) {
    locs <- str_locate_all(trans, "\\.|-") |>
        map_int(\(x) {
            if (nrow(x) == 0) return(NA_integer_)
            x[nrow(x),"end"]
        })
    out <- character(length(trans))
    out[is.na(locs)] <- trans[is.na(locs)]
    out[!is.na(locs)] <- str_sub(trans[!is.na(locs)], 1,
                                 locs[!is.na(locs)]-1)
    return(out)
}



# Lists to lookup HOGs for each gene and GO terms for HOG:
hog_go_list <- "~/_data/midgenomes/ortho-analyses/orthofinder-extraction/All_HOG_GO/N0-GO-by-HOG.tsv" |>
    read_tsv(col_types = "cc") |>
    mutate(go = str_split(go, ";")) |>
    (\(x) {
        z <- x$go
        names(z) <- x$hog
        return(z)
    })()
hog_gene_list <- "~/_data/midgenomes/ortho-analyses/orthofinder-extraction/All_HOG_GO/N0-GO-by-species-genes.tsv" |>
    read_tsv(col_types = "cccc") |>
    filter(species == "Tgraci") |>
    group_by(gene) |>
    summarize(hog = list(unique(hog))) |>
    (\(x) {
        z <- x$hog
        names(z) <- x$gene
        return(z)
    })()



gene_df <- "~/_data/Tgraci_genes.gff3" |>
    read_tsv(col_names = c("contig", "source", "type", "start", "end",
                           "score", "strand", "phase", "attributes"),
             comment = "#", col_types = "ccciicccc") |>
    mutate(contig = factor(contig, levels = paste0("CONTIG", 1:45)),
           gene = attributes |>
               str_remove_all("^ID|geneID|=") |>
               str_split(";") |>
               map_chr(\(x) x[[1]]),
           gene = trans2genes(gene)) |>
    select(contig, strand, start, end, gene) |>
    arrange(contig, strand, start)


# Functional annotation for all genes.
# (Have to use `read_lines` because it has different numbers of columns by row)
description_list <- "~/_data/midgenomes/go-terms/Tgraci_go-terms.tsv" |>
    read_lines(progress = FALSE) |>
    base::`[`(-1) |>
    str_split("\t") |>
    lapply(\(x) {
        descx <- x[startsWith(x, "description")]
        tibble(gene = x[[1]],
               desc = ifelse(length(descx) == 0, NA_character_,
                             tolower(paste(gsub("^description:", "", descx),
                                           collapse = ";"))))
    }) |>
    list_rbind() |>
    mutate(gene = trans2genes(gene)) |>
    filter(!is.na(desc)) |>
    # To remove duplicate descriptions:
    group_by(gene) |>
    summarize(desc = desc |>
                  str_split(";") |>
                  do.call(what = c) |>
                  unique() |>
                  list()) |>
    # convert to list for quick access:
    (\(x) {
        z <- x[["desc"]]
        names(z) <- x[["gene"]]
        return(z)
    })()



# GO terms for all Tgraci genes:
go_list <- "~/_data/midgenomes/go-terms/Tgraci_go-terms.tsv" |>
    read_lines(progress = FALSE) |>
    base::`[`(-1) |>
    str_split("\t") |>
    lapply(\(x) {
        gox <- x[startsWith(x, "go")]
        tibble(gene = x[[1]],
               go = ifelse(length(gox) == 0, NA_character_,
                           tolower(paste(gox, collapse = ";"))))
    }) |>
    list_rbind() |>
    mutate(gene = trans2genes(gene)) |>
    filter(!is.na(go)) |>
    # To remove duplicate GO terms and genes:
    group_by(gene) |>
    summarize(go = go |>
                  str_split(";") |>
                  do.call(what = c) |>
                  unique() |>
                  toupper() |>
                  list()) |>
    # Now add go terms from hogs:
    mutate(go = map2(gene, go, \(.gene, .go) {
        if (is.list(.go)) .go <- do.call(c, .go)
        .hog <- hog_gene_list[[.gene]]
        if (length(.hog) == 0) return(.go)
        hog_gos <- hog_go_list[[.hog]]
        .go <- unique(toupper(c(.go, hog_gos)))
        return(.go)
    })) |>
    # convert to list for quick access:
    (\(x) {
        z <- x[["go"]]
        names(z) <- x[["gene"]]
        return(z)
    })()




# Match contigs and positions to nearest gene it's within 5kb of.
# Arguments should be vectors, and it returns a list of dataframe to unnest
gene_matcher <- function(contig, pos, threshold = 5e3) {
    stopifnot(length(contig) == length(pos))
    map2(contig, pos, \(cg, ps) {
        cg_d <- gene_df[gene_df$contig == cg,]
        up <- (cg_d$strand == "+" & ps < cg_d$start &
                   ps > cg_d$start - threshold) |
            (cg_d$strand == "-" & ps > cg_d$end & ps < cg_d$end + threshold)
        down <- (cg_d$strand == "+" & ps > cg_d$end &
                     ps < cg_d$end + threshold) |
            (cg_d$strand == "-" & ps < cg_d$start & ps > cg_d$start - threshold)
        inside <- cg_d$start <= ps & cg_d$end >= ps
        stopifnot(!any((up & down) | (up & inside) | (down & inside)))
        lgl <- up | down | inside
        if (!any(lgl)) return(tibble(gene = NA, location = NA, dist = NA))
        ps_d <- cg_d[lgl,]
        ps_d[["location"]] <- ""
        ps_d[["location"]][up[lgl]] <- "upstream"
        ps_d[["location"]][down[lgl]] <- "downstream"
        ps_d[["location"]][inside[lgl]] <- "inside"
        ps_d[["dist"]] <- ifelse((ps_d$strand == "+" & ps_d$location == "upstream") |
                                     (ps_d$strand == "-" & ps_d$location == "downstream"),
                                 ps_d$start - ps, ps - ps_d$end)
        ps_d[["dist"]][inside[lgl]] <- -1
        return(ps_d[,c("gene", "location", "dist")])
    })
}









candi_genes_df <- candi_df |>
    mutate(gene_info = gene_matcher(contig, pos)) |>
    unnest(gene_info) |>
    arrange(desc(bf_db), mrk)


candi_genes_df



candi_funcs_df <- candi_genes_df |>
    mutate(desc = map(gene, \(g) {
        if (is.na(g)) return(NA_character_)
        if (! g %in% names(description_list)) return(NA_character_)
        description_list[[g]]
    }),
    go = map(gene, \(g) {
        if (is.na(g)) return(NA_character_)
        if (! g %in% names(go_list)) return(NA_character_)
        go_list[[g]]
    }))

candi_funcs_df |>
    select(contig, pos, bf_db, gene:go) |>
    mutate(across(desc:go, \(x) map_chr(x, \(xx) paste(xx, collapse = ";")))) |>
    select(-desc)




# =============================================================================*
# BLAST Results ----
# =============================================================================*

faa_file <- "~/Desktop/Tgraci_proteins.faa"

faa <- read.fasta(faa_file, "AA")
names(faa) <- names(faa) |> trans2genes()
for (i in 1:length(faa)) {
    attr(faa[[i]], "name") <- NULL
    attr(faa[[i]], "Annot") <- NULL
    attr(faa[[i]], "class") <- NULL
    faa[[i]] <- paste(faa[[i]], collapse = "")
}; rm(i)


candi_funcs_faa <- candi_funcs_df |>
    filter(!is.na(gene)) |>
    getElement("gene") |>
    unique() |>
    set_names() |>
    map_chr(\(x) faa[[x]])

# Now write to fasta file:
candi_funcs_faa |>
    imap(\(x, i) c(paste0(">", i), x)) |>
    unname() |>
    do.call(what = c) |>
    write_lines("~/Desktop/Tgraci-funct.faa")




#'
#' Installed blast stand-alone tools, then set `BLASTDB=~/blastdb`
#' inside `~/.bash_profile`
#'
#' Downloaded Drosophila melanogaster DB from Uniprot:
#' https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28proteome%3AUP000000803%29%29
#'
#' Then gunzipped it.
#'
#'
#' makeblastdb -in ~/Desktop/uniprotkb_proteome_UP000000803_2025_09_22.fasta \
#'     -dbtype prot -title Dmel_UniProt -taxid 7227 -out ${BLASTDB}/Dmel_UniProt
#'
#'
#' Then did search:
#'
#' blastp -query ~/Desktop/Tgraci-funct.faa \
#'      -db Dmel_UniProt \
#'      -out ~/Desktop/Tgraci-funct-Dmel.csv \
#'      -outfmt 20
#'
#'
#' Some of the top results:
#'
#'
#'

"~/Desktop/Tgraci-funct-Dmel.csv" |>
    read_csv(col_types = cols()) |>
    group_by(qaccver) |>
    filter(evalue == min(evalue)) |>
    ungroup() |>
    select(1:2, evalue) |>
    print(n = 30)


# # A tibble: 21 × 3
#    qaccver                saccver                   evalue
#    <chr>                  <chr>                      <dbl>
#  1 CONTIG18-_anno1.g4081  tr|Q9V3N4|Q9V3N4_DROME 4.1 e+  0
#  2 CONTIG37-_anno2.g6296  tr|Q8MLQ7|Q8MLQ7_DROME 1.67e-134
#  3 CONTIG37-_anno2.g6297  sp|Q9VZ64|6PGL_DROME   2.8 e+  0
#  4 CONTIG37+_anno1.g1347  tr|A8JNT4|A8JNT4_DROME 1.02e-170
#  5 CONTIG37+_anno1.g1348  tr|A8JNT6|A8JNT6_DROME 7.34e- 44
#  6 CONTIG44+_anno1.g3075  sp|Q9VL06|UFD4_DROME   0
#  7 CONTIG44+_anno1.g3078  tr|M9PFM5|M9PFM5_DROME 0
#  8 CONTIG44+_anno1.g3078  sp|P25439|BRM_DROME    0
#  9 CONTIG44+_anno1.g3078  tr|M9PFS6|M9PFS6_DROME 0
# 10 CONTIG44-_anno1.g3076  tr|M9PCA6|M9PCA6_DROME 0
# 11 CONTIG44-_anno1.g3076  tr|Q9VMJ5|Q9VMJ5_DROME 0
# 12 CONTIG44-_anno1.g3076  tr|Q9VGE7|Q9VGE7_DROME 0
# 13 CONTIG44-_anno2.g13176 tr|M9PCA6|M9PCA6_DROME 0
# 14 CONTIG44-_anno2.g13176 tr|Q9VMJ5|Q9VMJ5_DROME 0
# 15 CONTIG44-_anno2.g13176 tr|Q9VGE7|Q9VGE7_DROME 0
# 16 CONTIG45-_anno1.g9743  tr|Q9VH82|Q9VH82_DROME 2.03e- 21
# 17 CONTIG45-_anno1.g9744  tr|Q9VH81|Q9VH81_DROME 0
# 18 CONTIG45-_anno1.g9745  tr|O46050|O46050_DROME 0
# 19 CONTIG45-_anno1.g9746  tr|M9PC41|M9PC41_DROME 8.52e-  4
# 20 CONTIG45-_anno1.g9747  tr|A4V1P2|A4V1P2_DROME 2.1 e+  0
# 21 CONTIG45-_anno1.g9747  sp|Q9VSR3|ORB2_DROME   2.1 e+  0


#'
#' Most of these are pretty good matches, except for genes
#' `CONTIG18-_anno1.g4081`, `CONTIG37-_anno2.g6297`, and `CONTIG45-_anno1.g9747`
#'
#' Below are functions for each of the rest:
#'
#' CONTIG37-_anno2.g6296  Q8MLQ7
#'     - Enables D-glucose transmembrane transporter activity and trehalose
#'       transmembrane transporter activity.
#'       (https://www.alliancegenome.org/gene/FB:FBgn0034909)
#' CONTIG37+_anno1.g1347  A8JNT4
#'     - Regulates myosin phosphatase activity. Augments Ca2+ sensitivity of the contractile apparatus.
#'       (https://www.uniprot.org/uniprotkb/A8JNT4/entry)
#' CONTIG37+_anno1.g1348  A8JNT6
#'     - Regulates myosin phosphatase activity. Augments Ca2+ sensitivity of the contractile apparatus.
#'       (https://www.uniprot.org/uniprotkb/A8JNT4/entry)
#' CONTIG44+_anno1.g3075  Q9VL06
#'     -
#' CONTIG44+_anno1.g3078  M9PFM5
#'     -
#' CONTIG44+_anno1.g3078  P25439
#'     -
#' CONTIG44+_anno1.g3078  M9PFS6
#'     -
#' CONTIG44-_anno1.g3076  M9PCA6
#'     -
#' CONTIG44-_anno1.g3076  Q9VMJ5
#'     -
#' CONTIG44-_anno1.g3076  Q9VGE7
#'     -
#' CONTIG44-_anno2.g13176 M9PCA6
#'     -
#' CONTIG44-_anno2.g13176 Q9VMJ5
#'     -
#' CONTIG44-_anno2.g13176 Q9VGE7
#'     -
#' CONTIG45-_anno1.g9743  Q9VH82
#'     -
#' CONTIG45-_anno1.g9744  Q9VH81
#'     -
#' CONTIG45-_anno1.g9745  O46050
#'     -
#' CONTIG45-_anno1.g9746  M9PC41
#'     -

