#'
#' Check the 2026-11 BayPass results (`_data/baypass-qvalues.rds`, from step
#' 5 of `_R/baypass-pod-threshold.R`): per missing-data setting x
#' covariate, candidates at q <= 0.05 and q <= 0.10 (owner, 2026-11-07),
#' SNPs beyond the null (no POD SNP at least as extreme; p only a bound),
#' how consistently the strongest SNPs reach the BF ceiling across the 10
#' runs, and how much the strongest SNPs agree among settings and among
#' covariates.
#'
#' Agreement between two sets counts only SNPs tested in both (settings keep
#' different SNPs, and a SNP can fall in different subsamples): overlap of
#' the top 100 SNPs by empirical p (ties by BF, then z), and the Spearman
#' correlation of median BF.
#'
#' Output saved to `_data/baypass-check-candidates.txt`.
#'

library(tidyverse)

out_file <- "_data/baypass-check-candidates.txt"
sink(out_file, split = TRUE)

q_df <- read_rds("_data/baypass-qvalues.rds")
bf_ceiling <- 10 * log10((2 * 5000 - 1) * 0.99 / 0.01)



# =============================================================================*
# ---- per setting x covariate ----
# =============================================================================*


cat("== Per setting x covariate\n")
q_df |>
    group_by(tag, covariate) |>
    summarize(snps = n(),
              pi0 = round(pi0[[1]], 3),
              q05 = sum(q <= 0.05),
              q10 = sum(q <= 0.10),
              min_q = signif(min(q), 3),
              beyond_null = sum(beyond_null),
              median_at_ceiling = sum(bf_db >= bf_ceiling - 1e-3),
              .groups = "drop") |>
    as.data.frame() |>
    print()

cat("\n== SNPs beyond the null, by number of runs at the BF ceiling\n")
q_df |>
    filter(beyond_null) |>
    count(tag, covariate, n_runs_ceiling) |>
    pivot_wider(names_from = n_runs_ceiling, values_from = n, values_fill = 0,
                names_sort = TRUE) |>
    as.data.frame() |>
    print()

#' Neighboring SNPs share LD, so SNPs beyond the null are grouped into
#' regions: SNPs on the same contig within 1 kb of the previous one.
cat("\n== SNPs beyond the null, grouped into regions (SNPs within 1 kb)\n")
q_df |>
    filter(beyond_null) |>
    arrange(tag, covariate, contig, pos) |>
    group_by(tag, covariate, contig) |>
    mutate(region = cumsum(c(TRUE, diff(pos) > 1000))) |>
    group_by(tag, covariate, contig, region) |>
    summarize(from = min(pos), to = max(pos), snps = n(),
              max_z = round(max(z), 2), .groups = "drop") |>
    select(-region) |>
    as.data.frame() |>
    print()



# =============================================================================*
# ---- agreement ----
# =============================================================================*


agree <- function(x, y, k = 100) {
    #' x, y: rows of `q_df` for two sets
    both <- inner_join(x |> select(contig, pos, pa = p_emp, ba = bf_db, za = z),
                       y |> select(contig, pos, pb = p_emp, bb = bf_db, zb = z),
                       by = c("contig", "pos"))
    top_a <- both |> arrange(pa, desc(ba), desc(za)) |> slice_head(n = k)
    top_b <- both |> arrange(pb, desc(bb), desc(zb)) |> slice_head(n = k)
    tibble(tested_both = nrow(both),
           top100_shared = nrow(semi_join(top_a, top_b, by = c("contig", "pos"))),
           cor_bf = cor(both$ba, both$bb, method = "spearman"))
}

sets <- q_df |>
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
    mutate(ov = map2(a, b, \(i, j) agree(sets[[i]], sets[[j]]))) |>
    unnest(ov)

cat("\n== Agreement between settings (same covariate)\n")
pairs_df |>
    filter(cov_a == cov_b) |>
    select(covariate = cov_a, a = tag_a, b = tag_b, tested_both:cor_bf) |>
    mutate(cor_bf = round(cor_bf, 3)) |>
    as.data.frame() |>
    print()

cat("\n== Agreement between covariates (same setting)\n")
pairs_df |>
    filter(tag_a == tag_b) |>
    select(tag = tag_a, a = cov_a, b = cov_b, tested_both:cor_bf) |>
    mutate(cor_bf = round(cor_bf, 3)) |>
    as.data.frame() |>
    print()

sink()
#' Printed tables pad lines with spaces; strip them.
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
