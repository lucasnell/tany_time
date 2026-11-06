#'
#' Check the 2026-11 BayPass candidates (`_data/baypass-bf.rds`, from
#' `_R/baypass-analysis.R`): per missing-data setting x covariate, how
#' consistently each candidate SNP passes the 30 dB threshold across the 10
#' runs, how many sit at the BF ceiling, and how much candidate sets overlap
#' among settings and among covariates.
#'
#' Overlap between two sets counts only SNPs tested in both (settings keep
#' different SNPs, and a SNP can fall in different subsamples). Jaccard
#' index = shared candidates / candidates in either.
#'
#' Output saved to `_data/baypass-check-candidates.txt`.
#'

library(tidyverse)

out_file <- "_data/baypass-check-candidates.txt"
sink(out_file, split = TRUE)

bf_df <- read_rds("_data/baypass-bf.rds")
bf_ceiling <- 10 * log10((2 * 5000 - 1) * 0.99 / 0.01)



# =============================================================================*
# ---- per setting x covariate ----
# =============================================================================*


cat("== Candidates (median bf_db >= 30) per setting x covariate\n")
bf_df |>
    group_by(tag, covariate) |>
    summarize(snps = n(),
              runs_min = min(n_runs),
              runs_max = max(n_runs),
              max_bf_db = max(bf_db),
              candidates = sum(bf_db >= 30),
              cand_all_runs_30 = sum(bf_db >= 30 & n_runs_30 == n_runs),
              cand_median_at_ceiling = sum(bf_db >= bf_ceiling - 1e-3),
              .groups = "drop") |>
    as.data.frame() |>
    print()

cat("\n== Candidates by number of runs with bf_db >= 30\n")
bf_df |>
    filter(bf_db >= 30) |>
    count(tag, covariate, n_runs_30) |>
    pivot_wider(names_from = n_runs_30, values_from = n, values_fill = 0,
                names_sort = TRUE) |>
    as.data.frame() |>
    print()



# =============================================================================*
# ---- overlap ----
# =============================================================================*


overlap <- function(x, y) {
    #' x, y: rows of `bf_df` for two sets
    both <- inner_join(x |> select(contig, pos, a = bf_db),
                       y |> select(contig, pos, b = bf_db),
                       by = c("contig", "pos"))
    tibble(tested_both = nrow(both),
           cand_a = sum(both$a >= 30),
           cand_b = sum(both$b >= 30),
           shared = sum(both$a >= 30 & both$b >= 30),
           jaccard = shared / sum(both$a >= 30 | both$b >= 30),
           cor_bf = cor(both$a, both$b, method = "spearman"))
}

sets <- bf_df |>
    group_by(tag, covariate) |>
    group_split() |>
    (\(z) {
        set_names(z, map_chr(z, \(x) paste(x$tag[[1]], x$covariate[[1]])))
    })()

pairs_df <- combn(names(sets), 2) |>
    t() |>
    as_tibble(.name_repair = \(x) c("a", "b")) |>
    mutate(tag_a = word(a, 1), cov_a = word(a, 2),
           tag_b = word(b, 1), cov_b = word(b, 2)) |>
    filter(tag_a == tag_b | cov_a == cov_b) |>
    mutate(ov = map2(a, b, \(i, j) overlap(sets[[i]], sets[[j]]))) |>
    unnest(ov)

cat("\n== Overlap between settings (same covariate)\n")
cat("(cor_bf: Spearman correlation of median bf_db over SNPs tested in both)\n")
pairs_df |>
    filter(cov_a == cov_b) |>
    select(covariate = cov_a, a = tag_a, b = tag_b, tested_both:cor_bf) |>
    mutate(across(c(jaccard, cor_bf), \(x) round(x, 3))) |>
    as.data.frame() |>
    print()

cat("\n== Overlap between covariates (same setting)\n")
pairs_df |>
    filter(tag_a == tag_b) |>
    select(tag = tag_a, a = cov_a, b = cov_b, tested_both:cor_bf) |>
    mutate(across(c(jaccard, cor_bf), \(x) round(x, 3))) |>
    as.data.frame() |>
    print()

sink()
#' Printed tables pad lines with spaces; strip them.
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
