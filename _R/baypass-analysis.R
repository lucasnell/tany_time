

library(poolfstat)
library(seqinr)
library(clusterProfiler)
library(rrvgo)
library(treemap)
library(tidyverse)

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


beta_df <- crossing(b = 0L, r = 1:10) |>
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

# filt_blast_df |>
#     print(n = 100)





# # A tibble: 43 × 3
#    uniprot       evalue bf_db
#    <chr>          <dbl> <dbl>
#  1 M9PCM8     0          53.0
#  2 O46100     0          53.0
#  3 Q9VYK6     0          53.0
#  4 B4F7L9     6.06e-146  53.0
#  5 Q8MLQ7     1.67e-134  53.0
#  6 P48613     1.09e-123  53.0
#  7 Q9NI63     1.89e-118  53.0
#  8 A4V4D2     1.31e-113  53.0
#  9 P40792     2.29e-106  53.0
# 10 Q9VUD3     3.07e- 91  53.0
# 11 Q9VBU4     6.60e- 69  53.0
# 12 Q9VDD2     8.04e- 61  53.0
# 13 A0A0B4KGY5 6.21e- 20  53.0
# 14 P25724     1.80e- 12  53.0
# 15 Q9VFP2     5.10e- 11  53.0
# 16 B7Z0X5     7.20e-  8  53.0
# 17 M9PC41     2.23e-  6  53.0
# 18 Q7KUL1     2.51e-170  48.4
# 19 Q8MRC9     3.88e-101  45.5
# 20 Q9VQK4     3.22e- 49  45.5
# 21 D0Z756     0          42.9
# 22 Q9V4U7     2.04e-153  42.9
# 23 Q9V4U9     7.66e-126  42.9
# 24 P33270     1.99e- 11  42.9
# 25 A8JNT4     1.02e-170  39.7
# 26 A8JNT6     7.34e- 44  39.7
# 27 M9PGJ0     0          39.4
# 28 Q9VKD5     2.57e-106  39.4
# 29 A8DZ28     7.86e- 67  39.4
# 30 Q7JV61     1.07e- 37  39.4
# 31 P07207     3.63e- 50  37.9
# 32 Q9VPQ8     5.15e- 25  37.9
# 33 P17886     0          37.5
# 34 Q9VBV4     3.04e- 52  37.5
# 35 A0A0B4LHM9 5.96e- 44  37.5
# 36 Q9VQ47     0          37.0
# 37 B5RIV0     4.39e- 97  36.5
# 38 X2JKS9     3.51e- 50  36.5
# 39 Q9VHA5     1.72e- 26  36.5
# 40 Q9VBX4     6.40e- 18  36.5
# 41 Q9VK29     2.25e- 46  33.9
# 42 Q9VJA9     3.10e-155  32.4
# 43 A1ZAB3     4.6 e- 95  32.4


#'  1 M9PCM8     0          53.0
#'     - flight, larval visceral muscle development, and locomotion
#'     - https://www.uniprot.org/uniprotkb/M9PCM8/entry
#'  2 O46100     0          53.0
#'     - amino acid transmembrane transport, cell volume homeostasis,
#'       chloride ion homeostasis, chloride transmembrane transport,
#'       potassium ion homeostasis, potassium ion import across plasma membrane
#'     - https://www.uniprot.org/uniprotkb/O46100/entry
#'  3 Q9VYK6     0          53.0
#'     - long term (odor) memory
#'     - https://www.uniprot.org/uniprotkb/Q9VYK6/entry
#'  4 B4F7L9     6.06e-146  53.0
#'     - Chromosome mapping suggests that WDY is fully contained in the kl-1
#'       region, and WDY may correspond to this fertility factor.
#'     - https://www.uniprot.org/uniprotkb/B4F7L9/entry
#'  5 Q8MLQ7     1.67e-134  53.0
#'     - D-glucose transmembrane transport and trehalose transport
#'     - https://www.uniprot.org/uniprotkb/Q8MLQ7/entry
#'  6 P48613     1.09e-123  53.0
#'     - Enhances para sodium channel function. Required during pupal
#'       development to rescue adult paralysis and also protects adult flies
#'       against heat-induced lethality.
#'     - https://www.uniprot.org/uniprotkb/P48613/entry
#'  7 Q9NI63     1.89e-118  53.0
#'     - female and male meiotic nuclear division, spermatocyte division,
#'       spermatogonial cell division
#'     - https://www.uniprot.org/uniprotkb/Q9NI63/entry
#'  8 A4V4D2     1.31e-113  53.0
#'     - long-term (odor) memory
#'     - https://www.uniprot.org/uniprotkb/A4V4D2/entry
#'  9 P40792     2.29e-106  53.0
#'     - melanotic encapsulation of foreign target, memory,
#'       motor neuron axon guidance, muscle attachment,
#'       muscle cell development, myoblast fusion,
#'       myoblast proliferation, nephrocyte filtration,
#'       neuron projection development
#'     - https://www.uniprot.org/uniprotkb/P40792/entry
#' 10 Q9VUD3     3.07e- 91  53.0
#'     - brain development, neuron differentiation
#'     - https://www.uniprot.org/uniprotkb/Q9VUD3/entry





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


# ### NO LONGER USED
# # Create map to quickly pull out `goall` column for each `entrezid`:
# if (!file.exists("_data/goall_map.rds")) {
#     GO2ALLEGS <- as.list(org.Dm.eg.db::org.Dm.egGO2ALLEGS)
#     # Takes ~ 8 min:
#     goall_map <- keys(Dm_db) |>
#         set_names() |>
#         map(\(e) {
#             lgl <- map_lgl(GO2ALLEGS, \(x) e %in% x)
#             return(unique(names(GO2ALLEGS)[lgl]))
#         })
#     write_rds(goall_map, "_data/goall_map.rds")
#     rm(GO2ALLEGS)
# } else {
#     goall_map <- read_rds("_data/goall_map.rds")
# }




# columns(Dm_db)

all_prots <- db_select(Dm_db, keys = keys(Dm_db),
                       columns = c("ENTREZID", "UNIPROT", "GENENAME", "GO")) |>
    as_tibble() |>
    rename_with(tolower) |>
    filter(!is.na(uniprot)) |>
    # group_by(entrezid, uniprot, genename) |>
    # summarize(go = list(go), .groups = "drop") |>
    # mutate(goall = map(entrezid, \(e) goall_map[[e]]),
    #        go = map2(go, goall, \(x, y) {
    #            unique(c(x, y))
    #        })) |>
    # select(-goall, -entrezid) |>
    # unnest(go) |>
    filter(!is.na(go)) |>
    distinct(uniprot, go) |>
    mutate(term = suppressMessages(db_select(GO.db::GO.db, go, "TERM")[,"TERM"]))



term2gene <- all_prots |>
    select(go, uniprot)
term2name <- all_prots |>
    select(go, term)

genes <- filt_blast_df$uniprot[filt_blast_df$uniprot %in% all_prots$uniprot]


overrep <- enricher(genes,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    # minGSSize = 1,
                    # maxGSSize = 5,
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
    select(Description, p.adjust)

# # A tibble: 9 × 2
#   Description                                          p.adjust
#   <chr>                                                   <dbl>
# 1 regulation of locomotor rhythm                         0.0137
# 2 neuron projection morphogenesis                        0.0201
# 3 regulation of glycolytic process                       0.0201
# 4 behavioral response to starvation                      0.0251
# 5 Golgi to plasma membrane transport                     0.0302
# 6 germline ring canal formation                          0.0305
# 7 glial cell migration                                   0.0469
# 8 epidermis development                                  0.0469
# 9 regulation of compound eye photoreceptor development   0.0469


# Treemap figure (not necessary, probably)


or_bp_scores <- or_df |>
    filter(ont == "BP") |>
    mutate(geneID = geneID |> str_split("/")) |>
    unnest(geneID) |>
    mutate(bf_db = map_dbl(geneID, \(x) filt_blast_df$bf_db[filt_blast_df$uniprot == x][[1]])) |>
    group_by(ID) |>
    summarize(bf_db = max(bf_db)) |>
    (\(x) {z <- (x$bf_db); names(z) <- x$ID; return(z)})()

#' `org.Dm.eg.db` is the genome wide annotation for *Drosophila melanogaster*
sem_data <- GOSemSim::godata(annoDb = "org.Dm.eg.db",
                             ont = "BP", keytype = "ENTREZID")
or_bp_sim_mat <- calculateSimMatrix(names(or_bp_scores),
                                    semdata = sem_data,
                                    "org.Dm.eg.db", ont = "BP")

or_bp_red <- reduceSimMatrix(or_bp_sim_mat,
                             scores = or_bp_scores,
                             orgdb = "org.Dm.eg.db") |>
    as_tibble()
or_bp_red

treemap_p <- function() {
    .pal <- viridisLite::turbo(length(unique(or_bp_red$parent)), begin = 0.2)
    treemap(or_bp_red, index = c("parentTerm", "term"),
            vSize = "score", type = "index", title = "",
            lowerbound.cex.labels = 0.1,
            palette = .pal,
            fontcolor.labels = c("#FFFFFFDD", "#00000080"), bg.labels = 0,
            border.col = "#00000080")
}
treemap_p()


