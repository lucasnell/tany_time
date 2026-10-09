
#'
#' Write the files that HTCondor `.sub` files in `_bash/` queue jobs from.
#' Run from the repo root. Re-run if samples or QC decisions change.
#'
#' Outputs:
#'   _bash/samples/all-pools.txt            all sequenced pools (one per line)
#'   _bash/samples/all-pools-n-adults.txt   pool,n_adults (comma separated)
#'   _bash/samples/time-pools.txt           time-series pools passing QC
#'   _bash/samples/space-pools.txt          spatial pools passing QC
#'   _bash/03-select-time/baypass_beta_info.txt   ising_beta,run,seed
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
#' Pools passing QC (median coverage >= 30), copied from the output listed
#' in `_R/sample-qc.R` (section "Final samples to use").
#' Update these if the coverage cutoff changes (see TODO.md).
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

write_lines(time_pools, "_bash/samples/time-pools.txt")
write_lines(space_pools, "_bash/samples/space-pools.txt")


#'
#' BayPass covariate runs: Ising beta of 0 and 1, 10 runs each, with seeds.
#' This reproduces the seeds used for the 2025 runs
#' (originally written via `clipr::write_clip`; see `02-baypass-beta.sub`).
#'
set.seed(1323829402)
x <- expand.grid(beta = 0:1, run = 1:10)
x[["seed"]] <- sample.int(2^31-1, nrow(x))
write_csv(x, "_bash/03-select-time/baypass_beta_info.txt", col_names = FALSE)
