#'
#' Is diversity elevated around the CONTIG18 locus? (Owner, 2026-11-20;
#' exploratory.) Long-term balancing selection keeps linked neutral
#' variation old, so nucleotide diversity (π) is raised, and Tajima's D is
#' positive, in a region around the selected site, narrowing with
#' recombination. Per-pool npstat statistics in 10-kb windows
#' (`02-npstat`, `tany_dir("npstat")`; the 17 SN pools, and the 18
#' spatial pools of 2019, copied 2026-11-20) for the window
#' holding the locus (CONTIG18:370,001-380,000) and its neighbours, as
#' percentiles of each pool's genome-wide distribution (windows with
#' >= 5,000 unmasked sites).
#'
#' Caveat: npstat works from mapped reads, and reads carrying the B allele
#' map poorly (`_data/contig18-raw-assoc.txt`), so diversity at the locus
#' itself is underestimated; the neighbouring windows aren't affected.
#'
#' Writes `_data/contig18-diversity.txt`.
#'

library(tidyverse)

source("_R/data-dir.R")

out_file <- "_data/contig18-diversity.txt"
sink(out_file, split = TRUE)

ctg <- "CONTIG18"
w_loc <- 373124 %/% 10000 + 1
one_set <- function(set) {
    pools <- read_lines(sprintf("_bash/samples/%s-pools.txt", set))
    np <- map(pools, \(p) {
        data.table::fread(cmd = paste0("gzip -dc '",
                                       tany_dir("npstat", paste0(p, "_npstat_10kb.stat.gz")),
                                       "'"),
                          select = c("sequence", "window", "length", "Pi",
                                     "Tajima_D")) |>
            as_tibble() |>
            filter(length >= 5000) |>
            mutate(pool = p,
                   pi_pct = 100 * rank(Pi) / n(),
                   d_pct = 100 * rank(Tajima_D, na.last = "keep") /
                       sum(!is.na(Tajima_D)))
    }) |>
        list_rbind()
    near <- np |>
        filter(sequence == ctg, window >= w_loc - 5, window <= w_loc + 5)
    cat("\n== ", set, " pools (", length(pools), "): median over pools, ",
        ctg, " windows around the locus (window ", w_loc, " holds it)\n",
        sep = "")
    near |>
        group_by(window) |>
        summarize(start_kb = (window[[1]] - 1) * 10, pools = n(),
                  pi = signif(median(Pi), 3),
                  pi_percentile = round(median(pi_pct)),
                  tajima_d = round(median(Tajima_D, na.rm = TRUE), 2),
                  d_percentile = round(median(d_pct, na.rm = TRUE)),
                  .groups = "drop") |>
        as.data.frame() |>
        print()
    cat("Genome-wide medians over pools: pi", signif(median(np$Pi), 3),
        "; Tajima's D", round(median(np$Tajima_D, na.rm = TRUE), 2), "\n")
    cat("Locus window by pool:\n")
    near |>
        filter(window == w_loc) |>
        transmute(pool, pi = signif(Pi, 3), pi_percentile = round(pi_pct),
                  tajima_d = round(Tajima_D, 2), d_percentile = round(d_pct)) |>
        as.data.frame() |>
        print()
}
one_set("time")
one_set("space")

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
