#'
#' Are the PODs as dispersed as the real SNPs under the Ω the analysis
#' uses? (Owner, 2026-11-18, after the Ω correction test failed;
#' `_data/baypass-pod-correct.txt`.)
#'
#' The real covariate runs on subsample 1 and the POD covariate runs both
#' fix Ω at the subsample-1 Ω (`-omegafile`), so their XtX (BayPass's
#' measure of a SNP's differentiation across pools given Ω; posterior mean
#' `M_XtX`) are directly comparable. Only XtX is read from the real runs,
#' never Bayes factors.
#'
#' Criterion (fixed before reading any XtX): for each setting x covariate
#' (run 1 of each), the POD / real ratio of the median and of the 90th
#' percentile of M_XtX. Both within 0.95-1.05: the PODs match the data's
#' dispersion. Below 0.95: PODs are less dispersed, the null is too light,
#' and q-values would be optimistic. Above 1.05: the reverse. The upper
#' tail is not used, since selected SNPs can inflate the real one.
#' XtXst (BayPass's standardized XtX) is reported alongside for
#' information.
#'
#' Reads `tany_dir("baypass-beta")` and `tany_dir("baypass-pod")`. Writes
#' `_data/baypass-pod-xtx.txt`.
#'

library(tidyverse)

source("_R/data-dir.R")

out_file <- "_data/baypass-pod-xtx.txt"
sink(out_file, split = TRUE)

mp_tags <- c("mp050", "mp075", "mp100")
covariates <- c("hm_logN_6", "logN", "logN_lag1")

read_xtx <- function(tar, member) {
    tmp <- tempfile()
    untar(tar, files = member, exdir = tmp)
    data.table::fread(cmd = paste0("gzip -dc '", file.path(tmp, member), "'"),
                      select = c("M_XtX", "XtXst"), data.table = FALSE)
}

res <- crossing(tag = mp_tags, covariate = covariates) |>
    mutate(x = map2(tag, covariate, \(t, cv) {
        nr <- sprintf("baypass-%s-sub1-%s-beta0-run1", t, cv)
        real <- read_xtx(tany_dir("baypass-beta", paste0(nr, ".tar")),
                         file.path(nr, sprintf("tany-%s-sub1-%s-beta0-run1_summary_pi_xtx.out.gz",
                                               t, cv)))
        np <- sprintf("baypass-pod-%s-%s-run1", t, cv)
        pod <- read_xtx(tany_dir("baypass-pod", paste0(np, ".tar")),
                        file.path(np, sprintf("tany-pod-%s-%s-run1_summary_pi_xtx.out.gz",
                                              t, cv)))
        q <- \(x, p) unname(quantile(x, p))
        tibble(n_real = nrow(real), n_pod = nrow(pod),
               real_q50 = q(real$M_XtX, 0.5), pod_q50 = q(pod$M_XtX, 0.5),
               real_q90 = q(real$M_XtX, 0.9), pod_q90 = q(pod$M_XtX, 0.9),
               real_q99 = q(real$M_XtX, 0.99), pod_q99 = q(pod$M_XtX, 0.99),
               xtxst_real_q50 = q(real$XtXst, 0.5),
               xtxst_pod_q50 = q(pod$XtXst, 0.5))
    })) |>
    unnest(x) |>
    mutate(ratio_q50 = pod_q50 / real_q50,
           ratio_q90 = pod_q90 / real_q90,
           ratio_q99 = pod_q99 / real_q99,
           verdict = case_when(
               ratio_q50 >= 0.95 & ratio_q50 <= 1.05 &
                   ratio_q90 >= 0.95 & ratio_q90 <= 1.05 ~ "match",
               ratio_q50 < 0.95 | ratio_q90 < 0.95 ~ "POD less dispersed",
               TRUE ~ "POD more dispersed"))

cat("== M_XtX quantiles, real (subsample 1) vs. POD, run 1, same Omega\n")
res |>
    select(tag, covariate, n_real, n_pod, real_q50, pod_q50, real_q90,
           pod_q90, real_q99, pod_q99) |>
    mutate(across(where(is.double), \(x) signif(x, 4))) |>
    as.data.frame() |>
    print()
cat("\n== Ratios (POD / real) and verdict (median and 90th percentile",
    "within 0.95-1.05)\n")
res |>
    select(tag, covariate, ratio_q50, ratio_q90, ratio_q99, xtxst_real_q50,
           xtxst_pod_q50, verdict) |>
    mutate(across(where(is.double), \(x) round(x, 3))) |>
    as.data.frame() |>
    print()

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
