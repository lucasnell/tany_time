#'
#' Does the POD null hold for the real data? (Owner, 2026-11-19.)
#'
#' In BayPass's auxiliary covariate model, every SNP's BF depends on P(δ),
#' the share of associated SNPs estimated from all SNPs together; null
#' SNPs' BFs rise with it (`_data/baypass-pod-spike.txt`: P(δ) 0.021 with
#' 1% spiked vs. 0.0126 in the null POD; null median BF -3.8 vs. -6.4 dB).
#' So empirical p-values against the null POD are valid only if the real
#' data's P(δ) isn't above the POD's.
#'
#' Criterion (fixed before reading any real P(δ), owner's choice
#' 2026-11-19): per setting x covariate, the mean real P(δ) (posterior mean
#' M_P from `_summary_Pdelta.out`, over 3 subsamples x 10 runs) at or below
#' the highest P(δ) among the 10 POD runs: the POD null stands. Above it:
#' the null is too light for the real data, and no real BF is read until
#' the owner decides how to proceed.
#'
#' Only P(δ) is read from the real runs, never BFs.
#'
#' Reads `tany_dir("baypass-beta")` and `tany_dir("baypass-pod")`. Writes
#' `_data/baypass-pdelta-check.txt`.
#'

library(tidyverse)

source("_R/data-dir.R")

out_file <- "_data/baypass-pdelta-check.txt"
sink(out_file, split = TRUE)

mp_tags <- c("mp050", "mp075", "mp100")
covariates <- c("hm_logN_6", "logN", "logN_lag1")
n_runs <- 10

read_pdelta <- function(tar, member) {
    tmp <- tempfile()
    untar(tar, files = member, exdir = tmp)
    x <- read_table(file.path(tmp, member), col_types = cols())
    unlink(tmp, recursive = TRUE)
    x$M_P
}

real <- crossing(tag = mp_tags, sub = 1:3, covariate = covariates,
                 run = seq_len(n_runs)) |>
    mutate(p_delta = pmap_dbl(list(tag, sub, covariate, run), \(t, k, cv, r) {
        nm <- sprintf("baypass-%s-sub%i-%s-beta0-run%i", t, k, cv, r)
        read_pdelta(tany_dir("baypass-beta", paste0(nm, ".tar")),
                    file.path(nm, sprintf("tany-%s-sub%i-%s-beta0-run%i_summary_Pdelta.out",
                                          t, k, cv, r)))
    }))
pod <- crossing(tag = mp_tags, covariate = covariates, run = seq_len(n_runs)) |>
    mutate(p_delta = pmap_dbl(list(tag, covariate, run), \(t, cv, r) {
        nm <- sprintf("baypass-pod-%s-%s-run%i", t, cv, r)
        read_pdelta(tany_dir("baypass-pod", paste0(nm, ".tar")),
                    file.path(nm, sprintf("tany-pod-%s-%s-run%i_summary_Pdelta.out",
                                          t, cv, r)))
    }))

res <- real |>
    group_by(tag, covariate) |>
    summarize(real_mean = mean(p_delta), real_min = min(p_delta),
              real_max = max(p_delta), .groups = "drop") |>
    left_join(pod |>
                  group_by(tag, covariate) |>
                  summarize(pod_mean = mean(p_delta), pod_min = min(p_delta),
                            pod_max = max(p_delta), .groups = "drop"),
              by = c("tag", "covariate")) |>
    mutate(pod_null_holds = real_mean <= pod_max)

cat("== P(delta) (posterior mean), real (3 subsamples x 10 runs) vs. POD",
    "(10 runs)\n")
res |>
    mutate(across(where(is.double), \(x) signif(x, 3))) |>
    as.data.frame() |>
    print()
cat("\nBy subsample (real, mean over 10 runs):\n")
real |>
    group_by(tag, covariate, sub) |>
    summarize(p_delta = signif(mean(p_delta), 3), .groups = "drop") |>
    pivot_wider(names_from = sub, names_prefix = "sub", values_from = p_delta) |>
    as.data.frame() |>
    print()
cat("\nPOD null holds (real mean <= POD max) in", sum(res$pod_null_holds),
    "of", nrow(res), "setting x covariate pairs\n")

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
