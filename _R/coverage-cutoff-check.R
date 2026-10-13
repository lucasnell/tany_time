
#'
#' Per-pool coverage and read QC for the time-series (SN) pools, to inform
#' the median-coverage cutoff used to choose pools (`_R/sample-qc.R`).
#'
#' Inputs: `_data/bed-summs.rds` (from `_R/sample-qc.R`; coverage from
#' `01-variant-calling/04-bam-cov` BED files), MultiQC general stats, sample info, abundances.
#' Output: `_data/coverage-cutoff-check.csv`.
#' Assembly GC content is 32.5% (`_R/sample-qc.R`).
#'

library(tidyverse)

qc_dir <- "~/_data/_tany-time/trimmed/trimmed_reports/multiqc_data"

cov_df <- read_rds("_data/bed-summs.rds") |>
    # q20/q80: 20th/80th percentiles of per-site coverage; q20 >= 10 means
    # >= 80% of sites pass the per-site minimum coverage (10) in `01-snape`
    select(pool = sample, median_cov = median, mean_cov = mean,
           q20_cov = q20, q80_cov = q80, prop_zero_cov = p0)

# fastp/FastQC stats are per read file (R1, R2); average them per pool
mqc_df <- read_tsv(file.path(qc_dir, "multiqc_general_stats.txt"),
                   col_types = cols(), progress = FALSE) |>
    # fastp rows are named `<pool>_L002`; FastQC rows `trimmed_<pool>_L002_R*`
    mutate(pool = str_remove_all(Sample, "^trimmed_|_L002.*")) |>
    group_by(pool) |>
    summarize(pct_dup_fastp = mean(`fastp-pct_duplication`, na.rm = TRUE),
              pct_gc = mean(`fastqc-percent_gc`, na.rm = TRUE),
              .groups = "drop")

abund_df <- read_csv("_data/tany-abundances.csv", col_types = "icdd")

out_df <- read_csv("_data/full-sample-info.csv", col_types = "cDccidd") |>
    filter(str_detect(biotech_id, "^SN-[12]-")) |>
    rename(pool = biotech_id) |>
    mutate(info = str_remove(pool, "_S.*") |> str_split("-"),
           gen = map_chr(info, \(x) letters[as.integer(x[[2]])]),
           year = map_int(info, \(x) {
               y <- as.integer(x[[3]])
               ifelse(y < 70L, 2000L + y, 1900L + y)
           })) |>
    left_join(abund_df, by = c("year", "gen")) |>
    left_join(cov_df, by = "pool") |>
    left_join(mqc_df, by = "pool") |>
    mutate(n_haps = 2L * n_adults,
           passes_30 = median_cov >= 30) |>
    select(pool, year, gen, N, n_adults, n_haps, median_cov, mean_cov,
           q20_cov, q80_cov, prop_zero_cov, pct_dup_fastp, pct_gc, passes_30) |>
    arrange(median_cov)

stopifnot(!any(is.na(out_df$median_cov)), !any(is.na(out_df$pct_gc)))

write_csv(out_df, "_data/coverage-cutoff-check.csv")
print(out_df, n = Inf, width = Inf)
