
#'
#' Write the files that HTCondor `.sub` files in `_bash/` queue jobs from.
#' Run from the repo root. Re-run if samples or QC decisions change.
#'
#' Outputs:
#'   _bash/samples/all-pools.txt            all sequenced pools (one per line)
#'   _bash/samples/all-pools-n-adults.txt   pool,n_adults (comma separated)
#'   _bash/samples/time-pools-candidates.txt   all time-series (SN) pools
#'   _bash/samples/space-pools-candidates.txt  all spatial (2019) pools
#'   _bash/samples/time-pools-cov30.txt     time pools with median coverage
#'                                          >= 30 (the 2025 set; sensitivity)
#'   _bash/samples/space-pools-cov30.txt    spatial pools, same criterion
#'
#' The pools actually used (`time-pools.txt`, `space-pools.txt`) are written
#' by `_R/pool-qc.R` from QC metrics (`02-allele-freqs/03-pool-qc`).
#'   _bash/03-select-time/baypass_inputs_info.txt  tag,min_prop_pools
#'   _bash/03-select-time/baypass_init_info.txt    tag,subsample
#'   _bash/03-select-time/baypass_beta_info.txt    ising_beta,run,seed,tag
#'

library(tidyverse)


samp_df <- read_csv("_data/full-sample-info.csv", col_types = "cDccidd") |>
    arrange(biotech_id)

write_lines(samp_df$biotech_id, "_bash/samples/all-pools.txt")

samp_df |>
    mutate(x = paste(biotech_id, n_adults, sep = ",")) |>
    getElement("x") |>
    write_lines("_bash/samples/all-pools-n-adults.txt")


#'
#' Candidate pools: all time-series pools from site SN (generation 1 or 2)
#' and all spatial pools (collected in 2019).
#'
samp_df |>
    filter(str_detect(biotech_id, "^SN-[12]-")) |>
    getElement("biotech_id") |>
    write_lines("_bash/samples/time-pools-candidates.txt")
samp_df |>
    filter(format(date, "%Y") == "2019") |>
    getElement("biotech_id") |>
    write_lines("_bash/samples/space-pools-candidates.txt")


#'
#' Pools with median coverage >= 30, copied from the output listed in
#' `_R/sample-qc.R` (section "Final samples to use"). This was the 2025
#' inclusion rule; kept for sensitivity analyses.
#'
time_pools <- c("SN-1-02_S55", "SN-1-07_S23", "SN-1-08_S24", "SN-1-09_S25",
                "SN-1-10_S26", "SN-1-13_S56", "SN-1-14_S27", "SN-1-97_S61",
                "SN-1-99_S28", "SN-2-00_S29", "SN-2-02_S30", "SN-2-07_S31",
                "SN-2-08_S32", "SN-2-15_S34", "SN-2-92_S35", "SN-2-97_S37",
                "SN-2-99_S38")
space_pools <- c("Ash-19_S5", "Blik-19_S6", "Blo-19_S7", "Ellid-19_S8",
                 "Hrisatjorn_S11", "Lei-19_S13", "Lys-19_S14", "Mik-19_S15",
                 "MikW-19_S16", "MyBR-19_S17", "MyKS-19-B_S18", "MySN-19_S19",
                 "Ola-19_S20", "Pernuvatn_S21", "Skj-19_S22", "Sva-19_S39",
                 "Tjornin_S40", "Vik-19_S41")
stopifnot(all(c(time_pools, space_pools) %in% samp_df$biotech_id))

write_lines(time_pools, "_bash/samples/time-pools-cov30.txt")
write_lines(space_pools, "_bash/samples/space-pools-cov30.txt")


#'
#' BayPass is run for a grid of missing-data settings: minimum proportion
#' of pools with data at a SNP (decided 2026-10-13; all are reported, 0.75 is
#' the main analysis). Tags are used in all BayPass file names.
#'
mp_df <- tibble(min_prop_pools = c(0.5, 0.75, 1),
                tag = sprintf("mp%03i", round(100 * min_prop_pools)))
MAIN_TAG <- "mp075"
stopifnot(MAIN_TAG %in% mp_df$tag)
mp_df |>
    mutate(min_prop_pools = format(min_prop_pools, nsmall = 2)) |>
    select(tag, min_prop_pools) |>
    write_csv("_bash/03-select-time/baypass_inputs_info.txt",
              col_names = FALSE)

#' Omega runs: all SNPs (subsample 0) for every setting; thinned
#' subsamples 1-3 (Omega consistency check) for the main setting only.
bind_rows(tibble(tag = mp_df$tag, subsample = 0L),
          tibble(tag = MAIN_TAG, subsample = 1:3)) |>
    write_csv("_bash/03-select-time/baypass_init_info.txt", col_names = FALSE)

#' Covariate runs: Ising beta of 0 and 1, 10 runs each, for every setting.
#' Seeds reproduce those of the 2025 runs (originally written via
#' `clipr::write_clip`) and are the same across settings, so settings can
#' be compared run by run.
set.seed(1323829402)
x <- expand.grid(beta = 0:1, run = 1:10)
x[["seed"]] <- sample.int(2^31-1, nrow(x))
crossing(as_tibble(x), tag = mp_df$tag) |>
    arrange(tag, run, beta) |>
    write_csv("_bash/03-select-time/baypass_beta_info.txt", col_names = FALSE)
