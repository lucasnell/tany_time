
#'
#' Calibrate the pool-QC outlier rule (`_R/pool-qc.R`) on simulated clean
#' data, before any real QC metrics exist: how many pools does each rule
#' variant exclude when no pool is an outlier, and does it still catch
#' planted outliers?
#'
#' Clean data: multivariate normal with p = 5 metrics (moderate
#' correlations), n = 18 (space) or 30 (time) pools.
#' Planted: 2 pools shifted by 6 SD on 2 metrics (gross failure) and 1 pool
#' shifted by 4 SD on 1 metric (moderate).
#'
#' Output: `_data/pool-qc-rule-sim.csv`.
#'

library(tidyverse)

set.seed(2093148117)
n_reps <- 1000
p <- 5
Sigma <- matrix(0.3, p, p); diag(Sigma) <- 1

#' One rule variant applied to matrix X. Returns logical outlier vector.
#'   sd_only: use only score distance on all p PCs (robust Mahalanobis);
#'            otherwise SD on k PCs + OD (ROBPCA-style).
apply_rule <- function(X, q, sd_only, prop_var = 0.8) {
    Z <- apply(X, 2, \(x) (x - median(x)) / mad(x))
    mcd <- MASS::cov.rob(Z, method = "mcd")
    eig <- eigen(mcd$cov, symmetric = TRUE)
    k <- if (sd_only) p else
        min(which(cumsum(eig$values) / sum(eig$values) >= prop_var), p - 1L)
    Zc <- sweep(Z, 2, mcd$center)
    V <- eig$vectors[, 1:k, drop = FALSE]
    S <- Zc %*% V
    sd_dist <- sqrt(rowSums(sweep(S^2, 2, eig$values[1:k], "/")))
    out <- sd_dist > sqrt(qchisq(q, k))
    if (!sd_only) {
        od <- sqrt(rowSums((Zc - S %*% t(V))^2))
        od23 <- od^(2/3)
        out <- out | od > (median(od23) + mad(od23) * qnorm(q))^(3/2)
    }
    out
}

variants <- crossing(q = c(0.975, 0.99, 0.999), sd_only = c(FALSE, TRUE))

sim <- crossing(n = c(18L, 30L), rep = 1:n_reps) |>
    pmap(\(n, rep) {
        X <- MASS::mvrnorm(n, rep(0, p), Sigma)
        Xo <- X
        Xo[1:2, 1:2] <- Xo[1:2, 1:2] + 6
        Xo[3, 1] <- Xo[3, 1] + 4
        pmap(variants, \(q, sd_only) {
            clean <- apply_rule(X, q, sd_only)
            planted <- apply_rule(Xo, q, sd_only)
            tibble(n = n, rep = rep, q = q, sd_only = sd_only,
                   fp_clean = sum(clean),
                   gross_caught = sum(planted[1:2]),
                   moderate_caught = planted[3],
                   fp_planted = sum(planted[-(1:3)]))
        }) |>
            list_rbind()
    }, .progress = TRUE) |>
    list_rbind()

summ <- sim |>
    group_by(n, q, sd_only) |>
    summarize(mean_fp_clean = mean(fp_clean),
              prop_reps_any_fp = mean(fp_clean > 0),
              mean_gross_caught_of_2 = mean(gross_caught),
              prop_moderate_caught = mean(moderate_caught),
              mean_fp_with_outliers = mean(fp_planted),
              .groups = "drop")

write_csv(summ, "_data/pool-qc-rule-sim.csv")
print(summ, width = Inf)
