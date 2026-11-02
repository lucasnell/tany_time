#'
#' Check the revised time-series pool-QC rule (`_R/pool-qc.R`, 2026-11-02):
#' calibrated ROBPCA with `alpha = 0.5` instead of the default 0.75.
#' `alpha` is the share of pools ROBPCA fits its PCs to, so 0.75 tolerates
#' up to ~25% outliers and 0.5 up to ~50%. The first real run had 43% of
#' time pools far out on missing data, and the 0.75 rule kept them.
#'
#' Same designs as `_R/pool-qc-rule-sim.R` (p = 5 metrics, correlation 0.3,
#' bootstrap calibration to a 2.5% per-pool false-positive rate), with:
#'   - n = 21 (time pools left after the missing-data limit) and 30;
#'   - clean: no outliers;
#'   - planted: 2 pools 6 SD out on 2 metrics, 1 pool 4 SD out on 1 metric;
#'   - heavy: 30% of pools 6 SD out on 2 metrics.
#' Each rule's bootstrap uses its own `alpha`.
#'
#' Output: `_data/pool-qc-rule-sim-alpha.csv`.
#'

library(tidyverse)

set.seed(1460274813)
p <- 5
Sigma <- matrix(0.3, p, p); diag(Sigma) <- 1
n_reps <- 200
n_boot <- 300

robust_z <- function(X) apply(X, 2, \(x) (x - median(x)) / mad(x))
robpca_score <- function(Z, alpha) {
    pc <- rrcov::PcaHubert(Z, k = 0, scale = FALSE, alpha = alpha)
    od_ratio <- if (pc@cutoff.od > 0) pc@od / pc@cutoff.od else 0
    pmax(pc@sd / pc@cutoff.sd, od_ratio)
}
calibrated_rule <- function(X, alpha, fp_target = 0.025) {
    Z <- robust_z(X)
    R <- cov2cor(rrcov::getCov(rrcov::CovMcd(Z)))
    boot <- unlist(map(seq_len(n_boot), \(i) {
        robpca_score(robust_z(MASS::mvrnorm(nrow(Z), rep(0, ncol(Z)), R)),
                     alpha)
    }))
    robpca_score(Z, alpha) > quantile(boot, 1 - fp_target, names = FALSE)
}

grid <- crossing(n = c(21L, 30L), alpha = c(0.75, 0.5), rep = 1:n_reps) |>
    mutate(seed = sample.int(.Machine$integer.max, n()))

one_rep <- function(n, alpha, rep, seed) {
    set.seed(seed)
    X <- MASS::mvrnorm(n, rep(0, p), Sigma)
    Xp <- X
    Xp[1:2, 1:2] <- Xp[1:2, 1:2] + 6
    Xp[3, 1] <- Xp[3, 1] + 4
    n_heavy <- round(0.3 * n)
    Xh <- X
    Xh[1:n_heavy, 1:2] <- Xh[1:n_heavy, 1:2] + 6
    clean <- calibrated_rule(X, alpha)
    planted <- calibrated_rule(Xp, alpha)
    heavy <- calibrated_rule(Xh, alpha)
    tibble(n = n, alpha = alpha, rep = rep,
           fp_clean = sum(clean),
           gross_caught = sum(planted[1:2]),
           moderate_caught = planted[3],
           fp_planted = sum(planted[-(1:3)]),
           heavy_caught = sum(heavy[1:n_heavy]) / n_heavy,
           fp_heavy = sum(heavy[-(1:n_heavy)]))
}

sim <- parallel::mclapply(seq_len(nrow(grid)), \(i) {
    do.call(one_rep, as.list(grid[i, ]))
}, mc.cores = 8) |>
    list_rbind()

summ <- sim |>
    group_by(n, alpha) |>
    summarize(fp_rate_clean = mean(fp_clean) / n[1],
              mean_gross_caught_of_2 = mean(gross_caught),
              prop_moderate_caught = mean(moderate_caught),
              mean_fp_planted = mean(fp_planted),
              prop_heavy_caught = mean(heavy_caught),
              mean_fp_heavy = mean(fp_heavy),
              .groups = "drop")

write_csv(summ, "_data/pool-qc-rule-sim-alpha.csv")
print(summ, width = Inf)
