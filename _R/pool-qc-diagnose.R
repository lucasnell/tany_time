#'
#' Diagnose why the pool-QC rule (`_R/pool-qc.R`) keeps time-series pools
#' with nearly all sites masked (2026-11-02 run). Diagnostic only: it
#' doesn't change the rule or the pool lists.
#'
#' ROBPCA (`rrcov::PcaHubert`) fits its PCs to the h = alpha * n least
#' outlying pools (default alpha = 0.75), so it tolerates at most ~25%
#' outliers. This script counts how many pools are far from the bulk on
#' each metric, refits with alpha = 0.5 (the maximum ROBPCA allows), and
#' reports k, the h-subset, and the scores both ways. Scores are compared
#' with ROBPCA's own cutoffs (score > 1), not with the calibrated threshold.
#'
#' Reads `_data/pool-qc.csv` (from `_R/pool-qc.R`).
#' Output: `_data/pool-qc-diagnose.txt`.
#'

library(tidyverse)

metrics <- c("prop_missing", "n_private", "picard_dup", "fastp_gc",
             "target_rate")
# Same metrics per set as `_R/pool-qc.R`:
set_metrics <- list(time = metrics,
                    space = setdiff(metrics, "n_private"))

qc_df <- read_csv("_data/pool-qc.csv", show_col_types = FALSE)

robust_z <- function(X) {
    apply(X, 2, \(x) {
        s <- mad(x)
        if (s == 0) s <- sd(x)
        (x - median(x)) / s
    })
}

sink("_data/pool-qc-diagnose.txt", split = TRUE)

for (s in c("time", "space")) {
    m <- set_metrics[[s]]
    d <- qc_df |> filter(set == s, if_all(all_of(m), \(x) !is.na(x)))
    Z <- robust_z(as.matrix(d[, m]))
    cat(sprintf("\n==== %s: %i pools\n", s, nrow(d)))
    cat("Pools with |robust z| > 5 on each metric:\n")
    print(colSums(abs(Z) > 5))
    cat(sprintf("Pools with |robust z| > 5 on any metric: %i (%.0f%%)\n",
                sum(rowSums(abs(Z) > 5) > 0),
                100 * mean(rowSums(abs(Z) > 5) > 0)))
    fits <- map(c(0.75, 0.5), \(a) {
        set.seed(1702635893)
        pc <- rrcov::PcaHubert(Z, k = 0, scale = FALSE, alpha = a)
        od_ratio <- if (pc@cutoff.od > 0) pc@od / pc@cutoff.od else 0
        cat(sprintf("\nalpha = %.2f: k = %i, eigenvalues %s\n", a, pc@k,
                    paste(signif(pc@eigenvalues, 3), collapse = ", ")))
        tibble(pool = d$pool,
               !!sprintf("score_a%.2f", a) := pmax(pc@sd / pc@cutoff.sd,
                                                  od_ratio))
    })
    out <- d |>
        select(pool, prop_missing, n_private, target_rate,
               calibrated_include = include) |>
        left_join(fits[[1]], by = "pool") |>
        left_join(fits[[2]], by = "pool") |>
        arrange(desc(prop_missing))
    print(out, n = Inf, width = Inf)
}

sink()
