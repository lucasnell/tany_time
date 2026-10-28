#'
#' How many sites in the current time-series allele-frequency data are on
#' CONTIG34, the T. gracilentus assembly's own mitochondrial contig (it
#' carries COI MZ659759.1 at 99%; see `TODO.md`), and are any of them
#' BayPass candidates? Mitochondrial sites have far higher read depth than
#' nuclear ones and don't belong in a nuclear analysis.
#'
#' Reads `_data/allele-freqs-time.rds` and `_data/candi_df.rds` (both local
#' only). Writes `_data/mito-contig-sites.txt`.
#'


suppressPackageStartupMessages({
    library(tidyverse)
})

mt_contig <- "CONTIG34"

freqs <- readRDS("_data/allele-freqs-time.rds")
candi <- readRDS("_data/candi_df.rds")

by_type <- freqs |>
    mutate(type = ifelse(contig == mt_contig, mt_contig, "other contigs")) |>
    group_by(type) |>
    summarize(sites = n_distinct(contig, pos), pools = n_distinct(pool),
              median_reads = median(n_reads), .groups = "drop")

out <- c(paste("checked", Sys.Date()),
         paste("mitochondrial contig:", mt_contig),
         capture.output(print(as.data.frame(by_type), row.names = FALSE)),
         paste0("BayPass candidates on ", mt_contig, " (of ", nrow(candi),
                "): ", sum(candi$contig == mt_contig)))
writeLines(sub("[[:space:]]+$", "", out), "_data/mito-contig-sites.txt")
cat(out, sep = "\n")
