
#'
#' Per-pool mitochondrial contamination index from competitive mapping
#' (`_bash/01-variant-calling/00b-hologenome`, `02-bwa`).
#'
#' For each COI barcode group `bc_<code>` (other trap taxa; `bc_Procla` is
#' several Procladius species) the index is its read density relative to
#' the T. gracilentus mitogenome:
#'
#'   index = (reads on bc_<code> / bp of bc_<code>) /
#'           (reads on mt_Tgraci / bp of mt_Tgraci)
#'
#' Reads are MAPQ >= 20, primary. Mitochondrial copy number differs among
#' species and tissues, so the index is relative: compare it among pools
#' (e.g., against abundance), not as a fraction of DNA.
#'
#' Input: `species-counts/` from `pool-qc.tar`
#' (`_bash/02-allele-freqs/03-pool-qc`), extracted into `qc_dir`.
#' Output: `_data/mito-contam.csv`.
#'

library(tidyverse)

qc_dir <- "~/_data/_tany-time/pool-qc"

counts_df <- list.files(file.path(qc_dir, "species-counts"),
                        "_species-counts.tsv$", full.names = TRUE) |>
    set_names(\(x) str_remove(basename(x), "_species-counts.tsv$")) |>
    map(read_tsv, col_types = "cdd") |>
    list_rbind(names_to = "pool")

mt_df <- counts_df |>
    filter(group == "mt_Tgraci") |>
    transmute(pool, mt_reads = n_reads, mt_density = n_reads / length)

contam_df <- counts_df |>
    filter(str_starts(group, "bc_")) |>
    mutate(taxon = str_remove(group, "^bc_")) |>
    left_join(mt_df, by = "pool") |>
    mutate(density = n_reads / length,
           index = density / mt_density) |>
    select(pool, taxon, n_reads, length, mt_reads, index) |>
    arrange(taxon, desc(index))

stopifnot(!any(is.na(contam_df$mt_reads)))

write_csv(contam_df, "_data/mito-contam.csv")

contam_df |>
    select(pool, taxon, index) |>
    pivot_wider(names_from = taxon, values_from = index) |>
    print(n = Inf, width = Inf)
