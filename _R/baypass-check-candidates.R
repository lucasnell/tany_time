
#'
#' Check that `_data/candi_df.rds` matches the BayPass outputs currently in
#' `~/_data/_baypass/beta`, and how consistently each candidate SNP passes
#' the 30 dB threshold across the 10 independent runs.
#'
#' Output saved to `_data/baypass-check-candidates.txt`.
#'

library(tidyverse)

baypass_dir <- "~/_data/_baypass"

beta_df <- map(1:10, \(r) {
    sprintf("%s/beta/baypass-beta0-run%i/tany-beta0-run%i_summary_betai.out.gz",
            baypass_dir, r, r) |>
        read_table(col_types = "iidddd", skip = 1,
                   col_names = c("cov", "mrk", "m_beta", "sd_beta",
                                 "pip", "bf_db"),
                   progress = FALSE) |>
        mutate(run = r)
}) |>
    list_rbind()

med_df <- beta_df |>
    group_by(mrk) |>
    summarize(bf_db = median(bf_db), .groups = "drop")

candi_df <- read_rds("_data/candi_df.rds")
recomputed <- med_df$mrk[med_df$bf_db >= 30]

runs_df <- beta_df |>
    filter(mrk %in% recomputed) |>
    group_by(mrk) |>
    summarize(n_runs_30 = sum(bf_db >= 30), .groups = "drop")

out <- c(sprintf("SNPs tested: %i", nrow(med_df)),
         sprintf("Max bf_db in any run: %.2f", max(beta_df$bf_db)),
         sprintf("Candidates (median bf_db >= 30): %i", length(recomputed)),
         sprintf("Identical to _data/candi_df.rds: %s",
                 identical(sort(recomputed), sort(candi_df$mrk))),
         "Number of runs (of 10) with bf_db >= 30, per candidate:",
         capture.output(print(table(runs_df$n_runs_30))))

cat(out, sep = "\n")
write_lines(out, "_data/baypass-check-candidates.txt")
