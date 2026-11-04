#'
#' Write the BayPass covariate files for `_bash/03-select-time/03-baypass-beta`
#' (owner's choices, 2026-11-04). Each covariate runs in its own job,
#' because BayPass's AUX model with several covariates in one file fits a
#' multiple regression (manual v3.1, section 6.2), not one test per
#' covariate.
#'
#'   hm_logN_6  primary: log of the harmonic mean of (N + 1) over the 6
#'              generations (3 years) before the pool's generation. The
#'              long-run effective size of a fluctuating population is near
#'              the harmonic mean of its sizes; 6 generations is about the
#'              time N took to return above 10,000 after each crash.
#'   logN       sensitivity: log(1 + N) in the pool's generation (the
#'              covariate of the 2025 runs).
#'   logN_lag1  sensitivity: log(1 + N) in the generation before the pool's.
#'
#' N is the SN trap count (`_data/tany-abundances.csv`). Values are in the
#' order of `_bash/samples/time-pools.txt`, which is the pool order of the
#' BayPass inputs (`04-combine-syncs` names file); BayPass scales each
#' covariate to mean 0 and SD 1. Reads no allele frequencies.
#' Drafts and alternatives: `_R/crash-covariate-drafts.R`.
#'
#' Output: `_bash/03-select-time/covariates/<covariate>.txt` (one line,
#' space separated) and `_data/baypass-covariates.csv`.
#'

library(tidyverse)

abund <- read_csv("_data/tany-abundances.csv", show_col_types = FALSE) |>
    arrange(year, gen) |>
    mutate(g = row_number())

cov_df <- tibble(pool = read_lines("_bash/samples/time-pools.txt")) |>
    mutate(gen = if_else(str_sub(pool, 4, 4) == "1", "a", "b"),
           yy = as.integer(str_sub(pool, 6, 7)),
           year = if_else(yy < 70, 2000L + yy, 1900L + yy)) |>
    left_join(abund, by = c("year", "gen")) |>
    rowwise() |>
    mutate(hm_logN_6 = log(1 / mean(1 / (abund$N[abund$g >= g - 6 &
                                                     abund$g < g] + 1))),
           logN_lag1 = abund$logN[abund$g == g - 1]) |>
    ungroup() |>
    select(pool, year, gen, N, hm_logN_6, logN, logN_lag1)
stopifnot(!anyNA(cov_df), nrow(cov_df) == 18)

dir.create("_bash/03-select-time/covariates", showWarnings = FALSE)
for (v in c("hm_logN_6", "logN", "logN_lag1")) {
    cov_df[[v]] |>
        round(4) |>
        paste(collapse = " ") |>
        write_lines(sprintf("_bash/03-select-time/covariates/%s.txt", v))
}
write_csv(cov_df, "_data/baypass-covariates.csv")
print(mutate(cov_df, across(hm_logN_6:logN_lag1, \(x) round(x, 2))), n = Inf)
print(round(cor(select(cov_df, hm_logN_6:logN_lag1)), 2))
