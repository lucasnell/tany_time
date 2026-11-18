#'
#' Null-only comparison of two candidate test statistics for the BayPass
#' covariate runs (owner, 2026-11-18), before any real Bayes factor is read:
#'   BF   median over the 10 runs of BF (dB); capped at 59.96 dB (all 5,000
#'        retained samples include the SNP)
#'   z    median over the 10 runs of |M_Beta| / SD_Beta (posterior mean of
#'        the covariate effect over its posterior SD); no cap
#' plus the tie-break rule (BF, ties at the cap broken by z).
#'
#' On the PODs (no association by construction), per setting x covariate:
#'   1. resolution: null SNPs tied at the maximum; distinct values among
#'      the top 100;
#'   2. run-to-run stability: share of the variance among SNPs that is
#'      between SNPs rather than between runs (intraclass correlation,
#'      ICC), overall and among the top 1% by the median;
#'   3. agreement: Spearman correlation of the two medians; overlap of
#'      their top 100 and top 1,000 null SNPs.
#' Calibration isn't compared: empirical p-values against this null are
#' calibrated by construction for any statistic. Power needs SNPs with a
#' real association (`_R/baypass-pod-spike.R`).
#'
#' Reads `tany_dir("baypass-pod")`. Writes `_data/baypass-pod-stat-compare.txt`
#' and `tany_dir("baypass-pod", "pod-null-stats.rds")` (per POD SNP: median
#' BF and z, for `_R/baypass-pod-spike.R`).
#'

library(tidyverse)

source("_R/data-dir.R")

out_file <- "_data/baypass-pod-stat-compare.txt"
sink(out_file, split = TRUE)

mp_tags <- c("mp050", "mp075", "mp100")
covariates <- c("hm_logN_6", "logN", "logN_lag1")
n_runs <- 10
bf_ceiling <- 10 * log10((2 * 5000 - 1) * 0.99 / 0.01)

read_betai <- function(tag, cv, r) {
    nm <- sprintf("baypass-pod-%s-%s-run%i", tag, cv, r)
    member <- file.path(nm, sprintf("tany-pod-%s-%s-run%i_summary_betai.out.gz",
                                    tag, cv, r))
    tmp <- tempfile()
    untar(tany_dir("baypass-pod", paste0(nm, ".tar")), files = member,
          exdir = tmp)
    x <- data.table::fread(cmd = paste0("gzip -dc '", file.path(tmp, member),
                                        "'"),
                           select = c("M_Beta", "SD_Beta", "BF(dB)"))
    unlink(tmp, recursive = TRUE)
    #' SD_Beta = 0 when the SNP was never in the model in any sample
    #' (M_Beta = 0 too): z = 0, no evidence.
    z <- ifelse(x$SD_Beta > 0, abs(x$M_Beta) / x$SD_Beta, 0)
    n_sd0 <<- n_sd0 + sum(x$SD_Beta == 0)
    n_all <<- n_all + nrow(x)
    cbind(bf = x[["BF(dB)"]], z = z)
}
n_sd0 <- 0
n_all <- 0

#' Intraclass correlation: variance of SNP means over total variance
#' (one-way, n_runs per SNP).
icc <- function(m) {
    mu <- rowMeans(m)
    within <- mean(apply(m, 1, var))
    between <- var(mu) - within / ncol(m)
    between / (between + within)
}

null_stats <- list()
res <- list()
for (tag in mp_tags) {
    for (cv in covariates) {
        runs <- map(seq_len(n_runs), \(r) read_betai(tag, cv, r))
        bf <- sapply(runs, \(x) x[, "bf"])
        z <- sapply(runs, \(x) x[, "z"])
        rm(runs)
        med_bf <- apply(bf, 1, median)
        med_z <- apply(z, 1, median)
        null_stats[[paste(tag, cv)]] <- tibble(tag = tag, covariate = cv,
                                               bf = med_bf, z = med_z)
        top1_bf <- med_bf >= quantile(med_bf, 0.99)
        top1_z <- med_z >= quantile(med_z, 0.99)
        o_bf <- order(med_bf, med_z, decreasing = TRUE)
        o_z <- order(med_z, decreasing = TRUE)
        top_overlap <- function(k) length(intersect(o_bf[1:k], o_z[1:k])) / k
        res[[paste(tag, cv)]] <- tibble(
            tag = tag, covariate = cv,
            bf_at_ceiling = sum(med_bf >= bf_ceiling - 1e-6),
            bf_tied_at_max = sum(med_bf == max(med_bf)),
            z_tied_at_max = sum(med_z == max(med_z)),
            bf_distinct_top100 = n_distinct(sort(med_bf, decreasing = TRUE)[1:100]),
            z_distinct_top100 = n_distinct(sort(med_z, decreasing = TRUE)[1:100]),
            icc_bf = icc(bf), icc_z = icc(z),
            icc_bf_top1 = icc(bf[top1_bf, ]), icc_z_top1 = icc(z[top1_z, ]),
            spearman = cor(med_bf, med_z, method = "spearman"),
            overlap_top100 = top_overlap(100),
            overlap_top1000 = top_overlap(1000))
        rm(bf, z)
    }
}
res <- bind_rows(res)
cat("SNP x run results with SD_Beta = 0 (z set to 0):", n_sd0, "of", n_all,
    "\n\n")

cat("== 1. Resolution at the top of the null\n")
res |>
    select(tag, covariate, bf_at_ceiling, bf_tied_at_max, z_tied_at_max,
           bf_distinct_top100, z_distinct_top100) |>
    as.data.frame() |>
    print()
cat("\n== 2. Run-to-run stability (ICC: share of variance between SNPs)\n")
res |>
    select(tag, covariate, icc_bf, icc_z, icc_bf_top1, icc_z_top1) |>
    mutate(across(where(is.double), \(x) round(x, 3))) |>
    as.data.frame() |>
    print()
cat("\n== 3. Agreement between the medians (overlap: top k by BF with",
    "ties broken by z, vs. top k by z)\n")
res |>
    select(tag, covariate, spearman, overlap_top100, overlap_top1000) |>
    mutate(across(where(is.double), \(x) round(x, 3))) |>
    as.data.frame() |>
    print()

write_rds(bind_rows(null_stats), tany_dir("baypass-pod", "pod-null-stats.rds"),
          compress = "gz")

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
