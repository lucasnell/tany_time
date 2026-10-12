
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
#' Rule variants:
#'   mass_robpca:  base-R version (MASS::cov.rob MCD, then ROBPCA-style
#'                 score/orthogonal distance cutoffs); the first version of
#'                 `_R/pool-qc.R`
#'   mass_mahal:   MASS MCD robust Mahalanobis distance on all metrics
#'   rrcov_robpca: rrcov::PcaHubert (ROBPCA; Hubert et al. 2005) with its
#'                 built-in cutoffs; number of PCs chosen by PcaHubert
#'   rrcov_mahal:  rrcov::CovMcd (reweighted, consistency- and
#'                 small-sample-corrected MCD) robust Mahalanobis distance
#' each at cutoff quantiles q = 0.975, 0.99, 0.999.
#'
#' Outputs: `_data/pool-qc-rule-sim.csv` (fixed-cutoff variants) and
#' `_data/pool-qc-rule-sim-calibrated.csv` (calibrated rule).
#'

library(tidyverse)

set.seed(2093148117)
n_reps <- 1000
p <- 5
Sigma <- matrix(0.3, p, p); diag(Sigma) <- 1

#' One rule variant applied to matrix X. Returns logical outlier vector.
apply_rule <- function(X, q, method) {
    Z <- apply(X, 2, \(x) (x - median(x)) / mad(x))
    if (method == "rrcov_robpca") {
        pc <- rrcov::PcaHubert(Z, k = 0, scale = FALSE,
                               crit.pca.distances = q)
        return(!pc@flag)
    }
    if (method == "rrcov_mahal") {
        mcd <- rrcov::CovMcd(Z)
        d2 <- mahalanobis(Z, rrcov::getCenter(mcd), rrcov::getCov(mcd))
        return(d2 > qchisq(q, ncol(Z)))
    }
    sd_only <- method == "mass_mahal"
    prop_var <- 0.8
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

variants <- crossing(q = c(0.975, 0.99, 0.999),
                     method = c("mass_robpca", "mass_mahal", "rrcov_robpca",
                                "rrcov_mahal"))

sim <- crossing(n = c(18L, 30L), rep = 1:n_reps) |>
    pmap(\(n, rep) {
        X <- MASS::mvrnorm(n, rep(0, p), Sigma)
        Xo <- X
        Xo[1:2, 1:2] <- Xo[1:2, 1:2] + 6
        Xo[3, 1] <- Xo[3, 1] + 4
        pmap(variants, \(q, method) {
            clean <- apply_rule(X, q, method)
            planted <- apply_rule(Xo, q, method)
            tibble(n = n, rep = rep, q = q, method = method,
                   fp_clean = sum(clean),
                   gross_caught = sum(planted[1:2]),
                   moderate_caught = planted[3],
                   fp_planted = sum(planted[-(1:3)]))
        }) |>
            list_rbind()
    }, .progress = TRUE) |>
    list_rbind()

summ <- sim |>
    group_by(n, method, q) |>
    summarize(mean_fp_clean = mean(fp_clean),
              prop_reps_any_fp = mean(fp_clean > 0),
              mean_gross_caught_of_2 = mean(gross_caught),
              prop_moderate_caught = mean(moderate_caught),
              mean_fp_with_outliers = mean(fp_planted),
              .groups = "drop")

write_csv(summ, "_data/pool-qc-rule-sim.csv")
print(summ, width = Inf)



# =============================================================================*
# Calibrated ROBPCA (the rule used in `_R/pool-qc.R`) ----
# =============================================================================*
#'
#' Same clean/planted designs, fewer replicates because each one runs its
#' own bootstrap calibration (as `_R/pool-qc.R` does on real data).
#' Output: `_data/pool-qc-rule-sim-calibrated.csv`.
#'

robust_z <- function(X) apply(X, 2, \(x) (x - median(x)) / mad(x))
robpca_score <- function(Z) {
    pc <- rrcov::PcaHubert(Z, k = 0, scale = FALSE)
    od_ratio <- if (pc@cutoff.od > 0) pc@od / pc@cutoff.od else 0
    pmax(pc@sd / pc@cutoff.sd, od_ratio)
}
calibrated_rule <- function(X, fp_target = 0.025, n_boot = 300) {
    Z <- robust_z(X)
    R <- cov2cor(rrcov::getCov(rrcov::CovMcd(Z)))
    boot <- unlist(map(seq_len(n_boot), \(i) {
        robpca_score(robust_z(MASS::mvrnorm(nrow(Z), rep(0, ncol(Z)), R)))
    }))
    robpca_score(Z) > quantile(boot, 1 - fp_target, names = FALSE)
}

n_reps_cal <- 200
sim_cal <- crossing(n = c(18L, 30L), rep = 1:n_reps_cal) |>
    pmap(\(n, rep) {
        X <- MASS::mvrnorm(n, rep(0, p), Sigma)
        Xo <- X
        Xo[1:2, 1:2] <- Xo[1:2, 1:2] + 6
        Xo[3, 1] <- Xo[3, 1] + 4
        clean <- calibrated_rule(X)
        planted <- calibrated_rule(Xo)
        tibble(n = n, rep = rep, fp_clean = sum(clean),
               gross_caught = sum(planted[1:2]),
               moderate_caught = planted[3],
               fp_planted = sum(planted[-(1:3)]))
    }, .progress = TRUE) |>
    list_rbind()

summ_cal <- sim_cal |>
    group_by(n) |>
    summarize(fp_rate_clean = mean(fp_clean) / n[1],
              mean_fp_clean = mean(fp_clean),
              mean_gross_caught_of_2 = mean(gross_caught),
              prop_moderate_caught = mean(moderate_caught),
              mean_fp_with_outliers = mean(fp_planted),
              .groups = "drop")

write_csv(summ_cal, "_data/pool-qc-rule-sim-calibrated.csv")
print(summ_cal, width = Inf)
