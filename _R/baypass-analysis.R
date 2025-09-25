

library(poolfstat)
library(seqinr) # read.fasta
library(clusterProfiler) # enricher
library(tidyverse)

if (file.exists(".Rprofile")) source(".Rprofile")



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
    read_table(col_types = "cdcc", col_names = c("contig", "pos", "ref", "alt"),
               progress = FALSE) |>
    # Have to do it this way bc some are in scientific notation in file:
    mutate(pos = as.integer(pos)) |>
    select(contig, pos)


beta_df <- crossing(b = 0L, r = 1:10) |>
    pmap(\(b, r) {
        fn <- sprintf(paste0("%s/beta/baypass-beta%i-run%i/tany-beta%i-",
                             "run%i_summary_betai.out.gz"),
                      baypass_dir, b, r, b, r)
        fn |>
            read_table(col_types = "iidddd", skip = 1,
                       col_names = c("cov", "mrk", "m_beta", "sd_beta",
                                     "pip", "bf_db"),
                       progress = FALSE) |>
            mutate(run = r)

    }) |>
    list_rbind() |>
    select(run, everything())

bf_df <- beta_df |>
    group_by(mrk) |>
    summarize(bf_db = median(bf_db),
              bf = 10^(bf_db / 10)) |>
    (\(x) {
        # Checks to make sure it's okay to bind with `snp_df`:
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


# For use in other analyses:
if (!file.exists("_data/candi_df.rds")) write_rds(candi_df, "_data/candi_df.rds")







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
        out_idx <- which(ps_d[["dist"]] == min(ps_d[["dist"]]))[[1]]
        ps_d <- ps_d[out_idx,]
        return(ps_d[,c("gene", "location", "dist")])
    })
}









candi_genes_df <- candi_df |>
    mutate(gene_info = gene_matcher(contig, pos)) |>
    unnest(gene_info) |>
    arrange(desc(bf_db), mrk)


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
    mutate(across(desc:go, \(x) map_chr(x, \(xx) ifelse(isTRUE(is.na(xx)), xx,
                                                        paste(xx, collapse = ";"))))) |>
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

# # Now write to fasta file:
# candi_funcs_faa |>
#     imap(\(x, i) c(paste0(">", i), x)) |>
#     unname() |>
#     do.call(what = c) |>
#     write_lines("~/Desktop/Tgraci-funct.faa")




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

blast_df <- "~/Desktop/Tgraci-funct-Dmel.csv" |>
    read_csv(col_types = cols()) |>
    group_by(qaccver) |>
    filter(evalue == min(evalue)) |>
    filter(saccver == saccver[[1]]) |>
    ungroup() |>
    filter(evalue < 1e-4) |>
    mutate(uniprot = str_split(saccver, "\\|") |> map_chr(\(x) x[[2]]))

filt_blast_df <- blast_df |>
    select(1:2, evalue, uniprot) |>
    mutate(bf_db = map_dbl(qaccver, \(g) max(candi_funcs_df$bf_db[which(candi_funcs_df$gene == g)]))) |>
    group_by(uniprot) |>
    summarize(evalue = min(evalue), bf_db = max(bf_db), .groups = "drop") |>
    arrange(desc(bf_db), evalue)

filt_blast_df

# # A tibble: 8 × 3
#   uniprot    evalue bf_db
#   <chr>       <dbl> <dbl>
# 1 M9PCM8  0          53.0
# 2 Q9VYK6  0          53.0
# 3 B4F7L9  6.06e-146  53.0
# 4 Q9VUD3  3.07e- 91  53.0
# 5 Q9VBU4  6.60e- 69  53.0
# 6 Q9VDD2  8.04e- 61  53.0
# 7 P25724  1.80e- 12  53.0
# 8 Q9VFP2  5.10e- 11  53.0


#' 1 M9PCM8  0          53.0
#'     - flight, larval visceral muscle development, and locomotion
#'     - https://www.uniprot.org/uniprotkb/M9PCM8/entry
#' 2 Q9VYK6  0          53.0
#'     - long term (odor) memory
#'     - https://www.uniprot.org/uniprotkb/Q9VYK6/entry
#' 3 B4F7L9  6.06e-146  53.0
#'     - Chromosome mapping suggests that WDY is fully contained in the kl-1
#'       region, and WDY may correspond to this fertility factor.
#'     - https://www.uniprot.org/uniprotkb/B4F7L9/entry
#' 4 Q9VUD3  3.07e- 91  53.0
#'     - brain development, neuron differentiation
#'     - https://www.uniprot.org/uniprotkb/Q9VUD3/entry
#' 5 Q9VBU4  6.60e- 69  53.0
#'     - rRNA processing
#'     - https://www.uniprot.org/uniprotkb/Q9VBU4/entry
#' 6 Q9VDD2  8.04e- 61  53.0
#'     - behavioral response to starvation, cellular response to glucose
#'       starvation, cellular response to starvation, cholesterol homeostasis,
#'       lipid metabolic process, positive regulation of autophagy
#'     - https://www.uniprot.org/uniprotkb/Q9VDD2/entry
#' 7 P25724  1.80e- 12  53.0
#'     - dendrite morphogenesis, germ-line stem cell population maintenance,
#'       negative regulation of apoptotic process, negative regulation of
#'       synaptic assembly at neuromuscular junction, oogenesis,
#'       spermatogenesis
#'     - https://www.uniprot.org/uniprotkb/P25724/entry
#' 8 Q9VFP2  5.10e- 11  53.0
#'     - eye development, negative regulation of protein import into nucleus,
#'       positive regulation of apoptotic process, positive regulation of JNK
#'       cascade, positive regulation of protein catabolic process,
#'       proteasome-mediated ubiquitin-dependent protein catabolic process,
#'       protein destabilization, regulation of proteolysis
#'     - https://www.uniprot.org/uniprotkb/Q9VFP2/entry










# =============================================================================*
# =============================================================================*
# Over-representation test ----
# =============================================================================*
# =============================================================================*



# Doing this instead of loading `org.Dm.eg.db` bc that overrides dplyr::select
columns <- AnnotationDbi::columns
keys <- AnnotationDbi::keys
db_select <- AnnotationDbi::select
Dm_db <- org.Dm.eg.db::org.Dm.eg.db
Ontology <- BiocGenerics::Ontology





# columns(Dm_db)

all_prots <- db_select(Dm_db, keys = keys(Dm_db),
                       columns = c("ENTREZID", "UNIPROT", "GENENAME", "GO")) |>
    as_tibble() |>
    rename_with(tolower) |>
    filter(!is.na(uniprot), !is.na(go)) |>
    distinct(uniprot, go) |>
    mutate(term = suppressMessages(db_select(GO.db::GO.db, go, "TERM")[,"TERM"]))



term2gene <- all_prots |>
    select(go, uniprot)
term2name <- all_prots |>
    distinct(go, term)

genes <- filt_blast_df$uniprot[filt_blast_df$uniprot %in% all_prots$uniprot]


overrep <- enricher(genes,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    TERM2GENE=term2gene, TERM2NAME=term2name) |>
    as_tibble()
overrep


or_df <- overrep |>
    #' Ontology:
    #'   * MF - molecular function
    #'   * BP - biological process
    #'   * CC - cellular component
    mutate(ont = Ontology(ID))

or_df |>
    filter(ont == "BP") |>
    select(Description, p.adjust) |>
    print(n = 40)



# # A tibble: 35 × 2
#    Description                                        p.adjust
#    <chr>                                                 <dbl>
#  1 regulation of carbon utilization                     0.0241
#  2 negative regulation of protein import into nucleus   0.0241
#  3 positive regulation of gluconeogenesis               0.0241
#  4 endoplasmic reticulum localization                   0.0241
#  5 cellular response to glucose starvation              0.0241
#  6 protein destabilization                              0.0241
#  7 regulation of glycolytic process                     0.0241
#  8 triglyceride storage                                 0.0241
#  9 positive regulation of cell cycle                    0.0263
# 10 behavioral response to starvation                    0.0263
# 11 female germline ring canal stabilization             0.0263
# 12 mitochondrion localization                           0.0263
# 13 Golgi to plasma membrane transport                   0.0263
# 14 microtubule anchoring                                0.0263
# 15 protein K48-linked ubiquitination                    0.0280
# 16 regulation of proteolysis                            0.0280
# 17 nucleus localization                                 0.0280
# 18 nucleus organization                                 0.0280
# 19 larval visceral muscle development                   0.0280
# 20 cholesterol homeostasis                              0.0311
# 21 segmentation                                         0.0340
# 22 positive regulation of protein catabolic process     0.0340
# 23 flight                                               0.0354
# 24 anterior/posterior axis specification, embryo        0.0371
# 25 regulation of neuromuscular synaptic transmission    0.0375
# 26 rRNA processing                                      0.0417
# 27 nuclear migration                                    0.0417
# 28 eye development                                      0.0432
# 29 oocyte anterior/posterior axis specification         0.0432
# 30 ovarian nurse cell to oocyte transport               0.0433
# 31 autophagy                                            0.0452
# 32 positive regulation of JNK cascade                   0.0465
# 33 positive regulation of apoptotic process             0.0467
# 34 exocytosis                                           0.0500
# 35 somatic muscle development                           0.0500


