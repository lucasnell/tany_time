

library(tidyverse)
library(poolfstat)

source("_R/baypass-utils.R")
# These functions not necessary:
rm(compute_genetic_offset, compute.local.scores, fmd.dist, geno2YN, plot.omega, simulate.baypass, simulate.PCcorrelated.covariate)

#' Run the following to make file names and folders of the form
#' `*-run${r}-sub${i}` instead of `*-sub${i}-run${r}`
#' (this has since been corrected).
#' It assumes that runs are separated by folder (`run1`, `run2`, etc.).
#'
#' ```bash
#' for ((r=1; r<=3; r++)); do
#'     cd run${r}
#'     for ((i=1; i<=5; i++)); do
#'         cd baypass-beta_1.0-sub${i}-run${r}
#'         for f in tany-beta_1.0-sub${i}-run${r}*; do
#'             mv ${f} ${f/-sub${i}-run${r}/-run${r}-sub${i}}
#'         done
#'         cd ..
#'         mv baypass-beta_1.0-sub${i}-run${r} \
#'             baypass-beta_1.0-run${r}-sub${i}
#'     done
#'     cd ..
#' done
#'
#' # also move everything into a single folder:
#' for ((r=1; r<=3; r++)); do
#'     cd run${r}
#'     for ((i=1; i<=5; i++)); do
#'         cd baypass-beta_1.0-run${r}-sub${i}
#'         mv *.log *.out ../
#'         cd ..
#'         rm -r baypass-beta_1.0-run${r}-sub${i}
#'     done
#'     cd ..
#' done
#' ```
#'



bp_runs <- map(1:3, \(r) {
    concatenate_res(dir = sprintf("~/_data/_baypass/beta/run%i", r),
                    anaprefix = sprintf("tany-beta_1.0-run%i-sub", r),
                    anasep = "",
                    extension = "",
                    nsubsets = 5,
                    snpdet_prefix = "~/_data/_baypass/beta/tany.snpdet.sub",
                    retrieve_pi_xtx = FALSE,
                    retrieve_bfis = TRUE,
                    retrieve_c2 = FALSE) |>
        set_names(c("contig", "pos", "bf_db")) |>
        as_tibble() |>
        mutate(bf = 10^(bf_db / 10),
               run = factor(r, levels = 1:3)) |>
        select(run, everything())
}) |>
    list_rbind()


