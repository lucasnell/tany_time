#'
#' Read the 2026-11 BayPass covariate runs into one table of per-SNP
#' statistics (median over the 10 runs): Bayes factor (dB) and
#' z = |M_Beta| / SD_Beta (posterior mean of the covariate effect over its
#' posterior SD; z = 0 when SD_Beta = 0, i.e. the SNP was never in the
#' model). z breaks ties at the BF ceiling (owner, 2026-11-19; see
#' `_data/baypass-pod-stat-compare.txt` and `_data/baypass-pod-spike.txt`).
#'
#' Writes `_data/baypass-bf.rds`, read by `_R/baypass-pod-threshold.R`
#' (q-values), `_R/baypass-analysis.R`, and `_R/baypass-check-candidates.R`.
#'

library(tidyverse)

source("_R/data-dir.R")



# =============================================================================*
# Reading baypass output ----
# =============================================================================*

#'
#' 2026-11 rerun (`_bash/03-select-time/03-baypass-beta`): 3 missing-data
#' settings x 3 thinned subsamples x 3 covariates x 10 runs (Ising beta = 0),
#' one `.tar` per job in the local data folder (`tany_dir("baypass-beta")`),
#' each holding `_summary_betai.out.gz` and the subsample's SNP positions
#' (`tany.snpdet.sub<k>`). Subsamples are disjoint (every third SNP), so a
#' setting's SNPs are the 3 subsamples stacked. Per SNP, the Bayes factor
#' (dB) is the median over the 10 runs.
#'
#' Set `allow_partial <- TRUE` before sourcing to run on the jobs that have
#' finished (for drafting only; median over fewer runs).
#'

if (!exists("allow_partial")) allow_partial <- FALSE

beta_dir <- tany_dir("baypass-beta")
mp_tags <- c("mp050", "mp075", "mp100")
covariates <- c("hm_logN_6", "logN", "logN_lag1")
n_runs <- 10
#' Largest possible BF: PIP = 1 in all NVAL = 5,000 samples, with prior
#' P(delta = 1) mean 0.01 (BayPass reports 59.96 dB).
bf_ceiling <- 10 * log10((2 * 5000 - 1) * 0.99 / 0.01)

job_df <- crossing(tag = mp_tags, sub = 1:3, covariate = covariates,
                   run = 1:n_runs) |>
    mutate(name = sprintf("baypass-%s-sub%i-%s-beta0-run%i", tag, sub,
                          covariate, run),
           file = file.path(beta_dir, paste0(name, ".tar")),
           found = file.exists(file))

cat(sprintf("BayPass beta outputs found: %i of %i\n", sum(job_df$found),
            nrow(job_df)))
if (!all(job_df$found) && !allow_partial) {
    stop("Missing BayPass outputs; copy them to ", beta_dir, " or set ",
         "allow_partial <- TRUE to draft")
}

#' Extract one member of a job's tar to a temporary folder; returns its path.
untar_one <- function(name, member) {
    tmp <- tempfile()
    untar(file.path(beta_dir, paste0(name, ".tar")),
          files = file.path(name, member), exdir = tmp)
    file.path(tmp, name, member)
}

read_snpdet <- function(name, sub) {
    untar_one(name, paste0("tany.snpdet.sub", sub)) |>
        read_table(col_types = "cccc", progress = FALSE,
                   col_names = c("contig", "pos", "ref", "alt")) |>
        # Have to do it this way bc some are in scientific notation in file:
        mutate(pos = as.integer(as.numeric(pos)),
               mrk = row_number()) |>
        select(mrk, contig, pos)
}

read_betai <- function(name) {
    untar_one(name, sprintf("tany-%s_summary_betai.out.gz",
                            str_remove(name, "^baypass-"))) |>
        read_table(col_types = "iidddd", skip = 1,
                   col_names = c("cov", "mrk", "m_beta", "sd_beta",
                                 "pip", "bf_db"),
                   progress = FALSE) |>
        mutate(z = ifelse(sd_beta > 0, abs(m_beta) / sd_beta, 0))
}

#' Median over runs, per setting x subsample x covariate (read one group at
#' a time: all runs together would be ~75 million rows).
bf_df <- job_df |>
    filter(found) |>
    group_by(tag, sub, covariate) |>
    group_split() |>
    map(\(jd) {
        snps <- read_snpdet(jd$name[[1]], jd$sub[[1]])
        runs <- map(jd$name, read_betai) |>
            list_rbind()
        stopifnot(all(runs$mrk %in% snps$mrk))
        runs |>
            group_by(mrk) |>
            #' Counts over runs first: `summarize()` evaluates in order, so
            #' after `bf_db = median(bf_db)` later terms would see the
            #' median, not the runs.
            summarize(n_runs = n(),
                      n_runs_30 = sum(bf_db >= 30),
                      n_runs_ceiling = sum(bf_db >= bf_ceiling - 1e-3),
                      bf_db = median(bf_db),
                      z = median(z),
                      m_beta = median(m_beta),
                      .groups = "drop") |>
            (\(x) {
                # One row per SNP of the subsample:
                stopifnot(identical(x$mrk, snps$mrk))
                return(x)
            })() |>
            left_join(snps, by = "mrk") |>
            mutate(tag = jd$tag[[1]], sub = jd$sub[[1]],
                   covariate = jd$covariate[[1]])
    }) |>
    list_rbind() |>
    select(tag, covariate, sub, mrk, contig, pos, everything()) |>
    # To make sure they stay in order of the assembly:
    mutate(contig = factor(contig, levels = paste0("CONTIG", 1:45))) |>
    arrange(tag, covariate, contig, pos)

# For use in other analyses (all SNPs):
write_rds(bf_df, "_data/baypass-bf.rds", compress = "gz")
