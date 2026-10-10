
#'
#' Decide which pools to include from per-pool QC metrics
#' (`_bash/02-allele-freqs/03-pool-qc`), before any allele-frequency
#' analysis.
#'
#' Metrics: missing-data rate after masking, number of private SNPs,
#' duplication rate (Picard), GC content (fastp), and mapping rate.
#' This script deliberately doesn't read abundances or allele frequencies,
#' so inclusion can't depend on results.
#'
#' Steps:
#'   1. Download and extract `pool-qc.tar` into `qc_dir` (below).
#'   2. Run this script with `RULES_FINAL <- FALSE` to write the metric
#'      table and plots; look at them.
#'   3. Fill in `rules` (below), set `RULES_FINAL <- TRUE`, and re-run.
#'      This writes `_bash/samples/{time,space}-pools.txt`, used by
#'      `02-allele-freqs/04-combine-syncs`.
#'
#' Outputs: `_data/pool-qc.csv`, `_plots/pool-qc-{time,space}.pdf`, and
#' (once rules are final) the pool lists above.
#'

library(tidyverse)
library(patchwork)

qc_dir <- "~/_data/_tany-time/pool-qc"

#' Exclusion rules: a pool is excluded if any metric is outside
#' [min, max] (NA = no bound) or is missing. Set before looking at any
#' allele-frequency results, and record the rationale for each bound here.
RULES_FINAL <- FALSE
rules <- tribble(
    ~metric,        ~min, ~max,
    "prop_missing",   NA,   NA,
    "n_private",      NA,   NA,
    "picard_dup",     NA,   NA,
    "fastp_gc",       NA,   NA,
    "map_rate",       NA,   NA)



# =============================================================================*
# Read metrics ----
# =============================================================================*

samp_df <- read_csv("_data/full-sample-info.csv", col_types = "cDccidd") |>
    select(pool = biotech_id, n_adults) |>
    mutate(n_haps = 2L * n_adults)

files_df <- read_tsv(file.path(qc_dir, "pool-qc-files.tsv"),
                     col_types = "cddddd", na = "NA")

qc_df <- c("time", "space") |>
    map(\(s) {
        cands <- read_lines(sprintf("_bash/samples/%s-pools-candidates.txt",
                                    s))
        read_tsv(file.path(qc_dir, sprintf("pool-qc-sync_%s.tsv", s)),
                 col_types = "ciiddi") |>
            # Pools without SNAPE output have no row; keep them as NA
            right_join(tibble(pool = cands), by = "pool") |>
            mutate(set = s)
    }) |>
    list_rbind() |>
    left_join(files_df, by = "pool") |>
    left_join(samp_df, by = "pool") |>
    select(set, pool, n_haps, prop_missing, n_snps, n_private, picard_dup,
           fastp_dup, fastp_gc, map_rate, usable_rate)

stopifnot(identical(sort(rules$metric),
                    sort(c("prop_missing", "n_private", "picard_dup",
                           "fastp_gc", "map_rate"))))



# =============================================================================*
# Plots ----
# =============================================================================*

#' Metrics against each other, and a PCA of the scaled metrics (pools with
#' any missing metric are left out of the PCA).
plot_set <- function(s) {
    d <- qc_df |> filter(set == s)
    long <- d |>
        select(pool, all_of(rules$metric)) |>
        pivot_longer(-pool, names_to = "metric")
    p1 <- long |>
        ggplot(aes(value, reorder(pool, value))) +
        geom_point() +
        facet_wrap(~ metric, scales = "free_x", nrow = 1) +
        theme_classic() +
        theme(axis.title = element_blank(), axis.text.y = element_text(size = 6))
    d_pca <- d |>
        select(pool, all_of(rules$metric)) |>
        drop_na()
    pca <- prcomp(select(d_pca, -pool), scale. = TRUE)
    p2 <- as_tibble(pca$x[, 1:2]) |>
        mutate(pool = d_pca$pool) |>
        ggplot(aes(PC1, PC2, label = pool)) +
        geom_text(size = 2) +
        theme_classic() +
        labs(caption = sprintf("PC1 %.0f%%, PC2 %.0f%% of variance",
                               100 * summary(pca)$importance[2, 1],
                               100 * summary(pca)$importance[2, 2]))
    p1 / p2 + plot_annotation(title = paste("Pool QC:", s))
}

for (s in c("time", "space")) {
    ggsave(sprintf("_plots/pool-qc-%s.pdf", s), plot_set(s),
           width = 10, height = 9)
}



# =============================================================================*
# Apply rules ----
# =============================================================================*

fails <- function(x, lo, hi) {
    is.na(x) | (!is.na(lo) & x < lo) | (!is.na(hi) & x > hi)
}

reasons <- map(seq_len(nrow(rules)), \(i) {
    r <- rules[i, ]
    ifelse(fails(qc_df[[r$metric]], r$min, r$max), r$metric, NA_character_)
}) |>
    do.call(what = cbind) |>
    apply(1, \(x) paste(na.omit(x), collapse = ";"))

qc_df <- qc_df |>
    mutate(include = reasons == "",
           exclude_reasons = reasons)

write_csv(qc_df, "_data/pool-qc.csv")
print(qc_df, n = Inf, width = Inf)

if (!RULES_FINAL) {
    message("RULES_FINAL is FALSE: pool lists not written. ",
            "Set `rules`, then `RULES_FINAL <- TRUE`, and re-run.")
} else {
    if (all(is.na(rules$min) & is.na(rules$max))) {
        stop("RULES_FINAL is TRUE but no rule has a bound")
    }
    for (s in c("time", "space")) {
        qc_df |>
            filter(set == s, include) |>
            getElement("pool") |>
            write_lines(sprintf("_bash/samples/%s-pools.txt", s))
    }
    qc_df |>
        group_by(set) |>
        summarize(n_candidates = n(), n_included = sum(include)) |>
        print()
}
