#'
#' Bayes factor (BF) cutoff from pseudo-observed data (PODs;
#' `_R/baypass-pod-inputs.R` → `03-select-time/04-baypass-pod`). POD SNPs
#' have no association with any covariate, so their BFs from the covariate
#' model are a null distribution for the real runs (same Ω, pool sizes,
#' depths, missing data, filters, seeds, and MCMC settings).
#'
#' Reports, per missing-data setting x covariate:
#'   1. whether the POD mimics the data: FMD between the Ω re-estimated
#'      from the POD (core model) and the Ω it was simulated from, against
#'      FMD between real subsamples of a setting (0.22-0.33,
#'      `_data/baypass-omega-fmd.csv`);
#'   2. null quantiles of BF (dB) from single runs and of the median over
#'      the 10 runs (the statistic used for candidates);
#'   3. the null probability that the median BF reaches each of a grid of
#'      cutoffs, and the expected number of null SNPs passing it among the
#'      SNPs the real analysis tests (expected false positives).
#' Steps 1-3 use no real BF results. Comparing with the real candidate
#' counts (estimated false discovery rate) runs only with
#' `compare_real <- TRUE`, after the rule for choosing the cutoff is fixed
#' (`TODO.md`).
#'
#' Reads `tany_dir("baypass-pod")`. Writes `_data/baypass-pod-threshold.txt`
#' and `_data/baypass-pod-null.csv` (step 3 table).
#'

library(tidyverse)
library(geigen) # geigen function inside fmd.dist

source("_R/data-dir.R")

if (!exists("compare_real")) compare_real <- FALSE

out_file <- "_data/baypass-pod-threshold.txt"
sink(out_file, split = TRUE)

pod_dir <- tany_dir("baypass-pod")
mp_tags <- c("mp050", "mp075", "mp100")
covariates <- c("hm_logN_6", "logN", "logN_lag1")
n_runs <- 10
bf_ceiling <- 10 * log10((2 * 5000 - 1) * 0.99 / 0.01)
#' SNPs tested in the real analysis (all 3 subsamples; `00-baypass-inputs`
#' output, `_bash/RERUN-CHECKLIST.md` I1):
n_real <- c(mp050 = 1023857, mp075 = 955748, mp100 = 528372)
cutoffs <- c(10, 15, 20, 25, 30, 35, 40, 45, 50, bf_ceiling)

fmd.dist <- function(mat1, mat2) {
    return(sqrt(sum(log(geigen(mat1,mat2)$values)**2)))
}

#' Extract one member of a job's tar to a temporary folder; returns its path.
untar_one <- function(name, member) {
    tmp <- tempfile()
    untar(file.path(pod_dir, paste0(name, ".tar")),
          files = file.path(name, member), exdir = tmp)
    file.path(tmp, name, member)
}
read_omega <- function(fn) {
    read_table(fn, col_names = FALSE, col_types = cols()) |>
        as.matrix() |>
        unname()
}



# =============================================================================*
# ---- 1. does the POD mimic the data? ----
# =============================================================================*


cat("== 1. FMD: Omega re-estimated from POD vs. Omega simulated from\n")
tibble(tag = mp_tags) |>
    mutate(fmd = map_dbl(tag, \(t) {
        nm <- sprintf("baypass-pod-%s-core-run1", t)
        om_pod <- read_omega(untar_one(nm, sprintf("tany-pod-%s-core-run1_mat_omega.out", t)))
        om_sim <- read_omega(untar_one(nm, sprintf("tany-init-%s-sub1_mat_omega.out", t)))
        fmd.dist(om_pod, om_sim)
    })) |>
    as.data.frame() |>
    print()
cat("(Real subsamples within a setting: FMD 0.22-0.33.)\n")



# =============================================================================*
# ---- 2. null BF distribution ----
# =============================================================================*


job_df <- crossing(tag = mp_tags, covariate = covariates, run = 1:n_runs) |>
    mutate(name = sprintf("baypass-pod-%s-%s-run%i", tag, covariate, run),
           found = file.exists(file.path(pod_dir, paste0(name, ".tar"))))
cat(sprintf("\nPOD covariate-model outputs found: %i of %i\n",
            sum(job_df$found), nrow(job_df)))
stopifnot(all(job_df$found))

null_df <- job_df |>
    group_by(tag, covariate) |>
    group_split() |>
    map(\(jd) {
        map(jd$name, \(nm) {
            untar_one(nm, sprintf("%s_summary_betai.out.gz",
                                  str_replace(nm, "^baypass-pod", "tany-pod"))) |>
                read_table(col_types = "iidddd", skip = 1,
                           col_names = c("cov", "mrk", "m_beta", "sd_beta",
                                         "pip", "bf_db"),
                           progress = FALSE) |>
                select(mrk, bf_db)
        }) |>
            list_rbind() |>
            mutate(tag = jd$tag[[1]], covariate = jd$covariate[[1]])
    }) |>
    list_rbind()

med_df <- null_df |>
    group_by(tag, covariate, mrk) |>
    summarize(n_runs = n(), bf_db = median(bf_db), .groups = "drop")
stopifnot(all(med_df$n_runs == n_runs))

probs <- c(0.99, 0.999, 0.9999)
q_tbl <- function(x, statistic) {
    x |>
        group_by(tag, covariate) |>
        summarize(statistic = statistic,
                  snps = n(),
                  q99 = quantile(bf_db, probs[1]),
                  q999 = quantile(bf_db, probs[2]),
                  q9999 = quantile(bf_db, probs[3]),
                  max = max(bf_db),
                  at_ceiling = sum(bf_db >= bf_ceiling - 1e-3),
                  .groups = "drop")
}
cat("\n== 2. Null BF (dB) quantiles\n")
bind_rows(q_tbl(null_df, "single run"), q_tbl(med_df, "median of 10")) |>
    arrange(tag, covariate, statistic) |>
    mutate(across(q99:max, \(x) round(x, 2))) |>
    as.data.frame() |>
    print()



# =============================================================================*
# ---- 3. null exceedance and expected false positives ----
# =============================================================================*


#' p_null = proportion of POD SNPs with median BF >= cutoff (0 if none: an
#' upper bound is then ~3 / 200,000, the 95% bound for zero events);
#' expected_fp = p_null x SNPs tested in the real setting.
exceed_df <- crossing(med_df |> distinct(tag, covariate), cutoff = cutoffs) |>
    mutate(n_pod = map2_int(tag, covariate, \(t, cv) {
        sum(med_df$tag == t & med_df$covariate == cv)
    }),
    n_exceed = pmap_int(list(tag, covariate, cutoff), \(t, cv, ct) {
        sum(med_df$bf_db[med_df$tag == t & med_df$covariate == cv] >= ct - 1e-3)
    }),
    p_null = n_exceed / n_pod,
    expected_fp = p_null * n_real[tag])

write_csv(exceed_df, "_data/baypass-pod-null.csv")

cat("\n== 3. Median BF >= cutoff: null proportion and expected false positives\n")
cat("(SNPs tested: ", paste(names(n_real), n_real, sep = " = ", collapse = ", "),
    ")\n", sep = "")
exceed_df |>
    mutate(cutoff = round(cutoff, 2), p_null = signif(p_null, 3),
           expected_fp = round(expected_fp, 1)) |>
    select(tag, covariate, cutoff, n_exceed, p_null, expected_fp) |>
    as.data.frame() |>
    print()



# =============================================================================*
# ---- 4. comparison with real candidates (after the rule is fixed) ----
# =============================================================================*


#' Estimated false discovery rate at each cutoff: expected false positives
#' / real SNPs with median BF >= cutoff (`_data/baypass-bf.rds`).
if (compare_real) {
    bf_df <- read_rds("_data/baypass-bf.rds")
    stopifnot(all(bf_df$n_runs == n_runs))
    cat("\n== 4. Estimated FDR with real candidate counts\n")
    exceed_df |>
        mutate(n_real_pass = pmap_int(list(tag, covariate, cutoff), \(t, cv, ct) {
            sum(bf_df$bf_db[bf_df$tag == t & bf_df$covariate == cv] >= ct - 1e-3)
        }),
        fdr = expected_fp / pmax(n_real_pass, 1)) |>
        mutate(cutoff = round(cutoff, 2), expected_fp = round(expected_fp, 1),
               fdr = round(fdr, 3)) |>
        select(tag, covariate, cutoff, expected_fp, n_real_pass, fdr) |>
        as.data.frame() |>
        print()
}

sink()
#' Printed tables pad lines with spaces; strip them.
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
