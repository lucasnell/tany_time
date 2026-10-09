
#'
#' Checks of claims made in the October 2026 code reviews of the BayPass
#' time-series analysis. Output saved to `_data/review-checks.txt`.
#'
#'   A. Allele polarity of `freq` in `timeseries-data.csv` vs. the BayPass
#'      reference allele.
#'   B. Direction of association with logN for candidates, in both polarities.
#'   C. Genome-wide distribution of correlations with logN (random SNPs).
#'   D. Dependence of candidate correlations on the two 1997 pools.
#'   E. Number of SNPs >= 30 dB in each single run.
#'   F. SN pools sequenced vs. used, and abundance of unsampled crash years.
#'

library(tidyverse)

set.seed(1310488804)

baypass_dir <- "~/_data/_baypass"
out <- character(0)
say <- function(...) {
    x <- sprintf(...)
    cat(x, "\n")
    out <<- c(out, x)
}

pools <- read_lines("~/_data/_tany-time/snape/count_masked_time/count_masked_time.names.gz")
logN <- read_lines(paste0(baypass_dir, "/baypass-inputs/tany-log_n.txt")) |>
    str_split(" ") |>
    getElement(1) |>
    as.numeric()
stopifnot(length(pools) == length(logN))

# Reference-allele frequencies for all SNPs (rows = SNPs, cols = pools):
geno <- paste0(baypass_dir, "/baypass-inputs/tany.genobaypass") |>
    read_table(col_names = FALSE, col_types = cols(.default = "i"),
               progress = FALSE) |>
    as.matrix()
ref_cnt <- geno[, seq(1, ncol(geno), 2)]
tot_cnt <- ref_cnt + geno[, seq(2, ncol(geno), 2)]
ref_frq <- ref_cnt / tot_cnt
rm(geno, ref_cnt, tot_cnt)

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

candi_df <- read_rds("_data/candi_df.rds")
ts_df <- read_csv("_data/timeseries-data.csv", col_types = "icicicdiidd")
stopifnot(setequal(unique(ts_df$pool), pools))


# ---------------------------------------------------------------------------*
# A, B. Polarity and direction ----
# ---------------------------------------------------------------------------*

cand_df <- ts_df |>
    mutate(pool = factor(pool, levels = pools)) |>
    arrange(mrk, pool) |>
    group_by(mrk) |>
    summarize(r_freq_ref = cor(freq, ref_frq[mrk[[1]], ]),
              r_freq_logN = cor(freq, logN),
              r_ref_logN = cor(ref_frq[mrk[[1]], ], logN),
              .groups = "drop") |>
    left_join(beta_df |>
                  filter(mrk %in% candi_df$mrk) |>
                  group_by(mrk) |>
                  summarize(m_beta = median(m_beta), .groups = "drop"),
              by = "mrk")

say("A. Candidates: %i", nrow(cand_df))
say("A. freq tracks BayPass ref allele (r > 0.9): %i; alt allele (r < -0.9): %i; neither: %i",
    sum(cand_df$r_freq_ref > 0.9), sum(cand_df$r_freq_ref < -0.9),
    sum(abs(cand_df$r_freq_ref) <= 0.9))
say("B. freq vs logN negative: %i of %i", sum(cand_df$r_freq_logN < 0), nrow(cand_df))
say("B. ref-allele freq vs logN negative: %i of %i", sum(cand_df$r_ref_logN < 0), nrow(cand_df))
say("B. BayPass median m_beta negative: %i of %i", sum(cand_df$m_beta < 0), nrow(cand_df))
say("B. sign(m_beta) == sign(r_ref_logN) for all: %s",
    all(sign(cand_df$m_beta) == sign(cand_df$r_ref_logN)))

#' Polarize every SNP to the allele that is rarer on average across pools:
minor_frq <- function(f) if (mean(f) > 0.5) 1 - f else f
cand_minor_r <- map_dbl(cand_df$mrk, \(i) cor(minor_frq(ref_frq[i, ]), logN))
say("B. Mean-minor allele vs logN negative (candidates): %i of %i",
    sum(cand_minor_r < 0), length(cand_minor_r))


# ---------------------------------------------------------------------------*
# C. Genome-wide background ----
# ---------------------------------------------------------------------------*

rand_i <- sample.int(nrow(ref_frq), 1e5)
rand_r_ref <- apply(ref_frq[rand_i, ], 1, cor, y = logN)
rand_r_minor <- apply(ref_frq[rand_i, ], 1, \(f) cor(minor_frq(f), logN))
# Expected P(|r| > 0.7) for n = 17 under independence:
n <- length(logN)
t07 <- 0.7 * sqrt((n - 2) / (1 - 0.7^2))
say("C. Random SNPs: %i", length(rand_i))
say("C. P(|r| > 0.7), observed: %.4f; expected under independence: %.4f",
    mean(abs(rand_r_ref) > 0.7), 2 * pt(-t07, n - 2))
say("C. Mean-minor allele vs logN negative (random SNPs): %.3f; mean r = %.3f",
    mean(rand_r_minor < 0), mean(rand_r_minor))
say("C. Binomial P(>= %i of %i negative | p = %.3f): %.2e",
    sum(cand_minor_r < 0), length(cand_minor_r), mean(rand_r_minor < 0),
    pbinom(sum(cand_minor_r < 0) - 1, length(cand_minor_r),
           mean(rand_r_minor < 0), lower.tail = FALSE))


# ---------------------------------------------------------------------------*
# D. 1997 pools ----
# ---------------------------------------------------------------------------*

p97 <- which(str_detect(pools, "-97_"))
say("D. 1997 pools: %s; logN = %s", paste(pools[p97], collapse = ", "),
    paste(round(logN[p97], 2), collapse = ", "))
say("D. logN rank of 1997 pools (1 = lowest): %s",
    paste(rank(logN)[p97], collapse = ", "))
abs_r_all <- map_dbl(cand_df$mrk, \(i) abs(cor(ref_frq[i, ], logN)))
abs_r_no97 <- map_dbl(cand_df$mrk, \(i) abs(cor(ref_frq[i, -p97], logN[-p97])))
say("D. Candidates, median |Pearson r| with logN: all pools %.2f; without 1997 %.2f",
    median(abs_r_all), median(abs_r_no97))
r_crit15 <- qt(0.975, 13) / sqrt(13 + qt(0.975, 13)^2)
say("D. Candidates with |r| > %.2f (nominal P < 0.05, n = 15) without 1997: %i of %i",
    r_crit15, sum(abs_r_no97 > r_crit15), length(abs_r_no97))


# ---------------------------------------------------------------------------*
# E. Single-run counts ----
# ---------------------------------------------------------------------------*

run_counts <- beta_df |>
    group_by(run) |>
    summarize(n30 = sum(bf_db >= 30), .groups = "drop")
say("E. SNPs >= 30 dB per run: %s", paste(run_counts$n30, collapse = ", "))
say("E. SNPs >= 30 dB in at least one run: %i",
    n_distinct(beta_df$mrk[beta_df$bf_db >= 30]))


# ---------------------------------------------------------------------------*
# F. Sampling coverage of crashes ----
# ---------------------------------------------------------------------------*

abund_df <- read_csv("_data/tany-abundances.csv", col_types = "icdd")
sn_df <- read_csv("_data/full-sample-info.csv", col_types = "cDccidd") |>
    # Time-series pools only (excludes 2019 spatial pool `MySN-19`):
    filter(str_detect(biotech_id, "^SN-[12]-")) |>
    mutate(info = str_remove(biotech_id, "_S.*") |> str_split("-"),
           gen = map_chr(info, \(x) letters[as.integer(x[[2]])]),
           year = map_int(info, \(x) {
               y <- as.integer(x[[3]])
               ifelse(y < 70L, 2000L + y, 1900L + y)
           }),
           used = biotech_id %in% pools) |>
    left_join(abund_df, by = c("year", "gen"))
say("F. SN pools sequenced: %i; used: %i", nrow(sn_df), sum(sn_df$used))
say("F. Not used: %s", sn_df |> filter(!used) |> arrange(year, gen) |>
        mutate(x = sprintf("%i%s (N=%.0f)", year, gen, N)) |>
        getElement("x") |> paste(collapse = ", "))
say("F. Min N among used pools: %.0f", min(sn_df$N[sn_df$used]))
say("F. Generations 1977-2015 with N < 100: %s", abund_df |>
        filter(N < 100) |>
        mutate(x = sprintf("%i%s (%.0f)", year, gen, N)) |>
        getElement("x") |> paste(collapse = ", "))

write_lines(out, "_data/review-checks.txt")
