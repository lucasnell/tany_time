
#'
#' Decide which pools to include from per-pool QC metrics
#' (`_bash/02-allele-freqs/03-pool-qc`), before any allele-frequency
#' analysis.
#'
#' Metrics: missing-data rate after masking, number of private SNPs,
#' duplication rate (Picard), GC content (fastp), and mapping rate to
#' T. gracilentus under competitive mapping.
#' This script deliberately doesn't read abundances or allele frequencies,
#' so inclusion can't depend on results.
#'
#' Rule for the time series (decided 2026-10-10, before any QC metrics
#' existed; revised 2026-11-02, see below): a pool is excluded if any metric
#' is missing, if more than half its sites are masked, or if it is an outlier
#' in a robust PCA of the QC metrics.
#'
#'   0. Exclude pools with `prop_missing` > MAX_MISSING (0.5): they lack
#'      data at most sites.
#'   1. Standardize each metric robustly over the remaining pools:
#'      (x - median) / MAD.
#'   2. ROBPCA (`rrcov::PcaHubert`; Hubert et al. 2005,
#'      doi:10.1198/004017004000000563) fitted to the ALPHA = 0.5 least
#'      outlying pools, with its default number of PCs and default (0.975)
#'      cutoffs for score distance (SD; within the PCs) and orthogonal
#'      distance (OD; off the PCs).
#'   3. Outlier score = max(SD / SD cutoff, OD / OD cutoff).
#'   4. Calibrate: ROBPCA's cutoffs flag 4-9% of pools in clean simulated
#'      data at n = 18-30 (`_R/pool-qc-rule-sim.R`), so the threshold for
#'      the score is set by parametric bootstrap. Simulate clean data
#'      (multivariate normal, same n, robust correlation of the metrics from
#'      reweighted MCD, `rrcov::CovMcd`), apply steps 1-3, and use the
#'      (1 - FP_TARGET) quantile of simulated scores. This targets a
#'      per-pool false-positive rate of FP_TARGET under the clean model.
#'   5. Outlier if score > threshold.
#'
#' Revision (owner, 2026-11-02): the first run (no step 0, ROBPCA's default
#' alpha = 0.75) kept time pools with every site masked. 13 of 30 time pools
#' were far out on missing data, more than the ~25% outliers alpha = 0.75
#' tolerates (`_R/pool-qc-diagnose.R`). Steps 0 and alpha = 0.5 were added
#' (`_R/pool-qc-rule-sim-alpha.R` checks the calibrated rule at alpha =
#' 0.5). `_data/pool-qc-sensitivity.csv` gives the time pools kept by the
#' first-run rule and by other missing-data limits.
#'
#' Override (owner, 2026-11-02): SN-1-95 is included although the rule
#' excludes it. Only its `prop_missing` is off (0.237 vs. 0.163-0.175 in the
#' kept pools; its other metrics are within their range), and it passes the
#' 0.5 limit; it scores high because `prop_missing` varies so little among
#' good pools. No general rule change brought it back without also keeping
#' clearly poor pools. The rule's own list (without SN-1-95) is reported as
#' a sensitivity analysis (`include_rule` in `_data/pool-qc.csv`).
#'
#' Spatial pools get no QC (owner, 2026-11-02): they were processed right
#' after collection, not stored for years with other species, so they
#' shouldn't be degraded. All candidates are included. Their metrics are
#' still saved.#'
#' Steps:
#'   1. Download and extract `pool-qc.tar` into `qc_dir` (below).
#'   2. Run this script. It writes `_bash/samples/{time,space}-pools.txt`,
#'      used by `02-allele-freqs/04-combine-syncs`.
#'
#' Outputs: `_data/pool-qc.csv`, `_data/pool-qc-sensitivity.csv`,
#' `_plots/pool-qc-time.pdf`, and the pool lists above.
#'

library(tidyverse)
library(patchwork)

source("_R/data-dir.R")
qc_dir <- tany_dir("pool-qc")

#' QC metrics used in the rule (definitions in `03-pool-qc.sh` and
#' `pool-qc-files.py`). The mapping rate is `target_rate`: reads assigned to
#' T. gracilentus contigs (MAPQ >= 20, primary) under competitive mapping,
#' over all reads, excluding reads on a mitochondrial-like stretch of
#' CONTIG37 (`mitolike_rate`; decided 2026-10-30, before any QC metrics
#' existed, because they are mitochondrial reads).
metrics <- c("prop_missing", "n_private", "picard_dup", "fastp_gc",
             "target_rate")

#' Rule parameters (fixed before seeing any metrics, except MAX_MISSING and
#' ALPHA, set 2026-11-02; see above):
FP_TARGET <- 0.025     # per-pool false-positive rate under the clean model
N_BOOT <- 2000         # simulated data sets for calibration
RULE_SEED <- 1702635893 # MCD and bootstrap use random numbers
MAX_MISSING <- 0.5     # step 0
ALPHA <- 0.5           # share of pools ROBPCA fits its PCs to
OVERRIDE_INCLUDE <- c("SN-1-95_S59")  # owner's override; see above



# =============================================================================*
# Read metrics ----
# =============================================================================*

samp_df <- read_csv("_data/full-sample-info.csv", col_types = "cDccidd") |>
    select(pool = biotech_id, n_adults) |>
    mutate(n_haps = 2L * n_adults)

files_df <- read_tsv(file.path(qc_dir, "pool-qc-files.tsv"),
                     col_types = "cdddddddd", na = "NA")

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
           fastp_dup, fastp_gc, target_rate, nontarget_rate, map_rate,
           usable_rate, mitolike_rate)




# =============================================================================*
# Robust PCA outlier rule ----
# =============================================================================*

#' Robust standardization of each column: (x - median) / MAD.
robust_z <- function(X) {
    apply(X, 2, \(x) {
        s <- mad(x)
        if (s == 0) s <- sd(x)
        if (s == 0) stop("metric with no variation")
        (x - median(x)) / s
    })
}

#' ROBPCA outlier scores: max(SD / SD cutoff, OD / OD cutoff), plus k and
#' the first two robust PC scores.
robpca_scores <- function(Z, alpha) {
    pc <- rrcov::PcaHubert(Z, k = 0, scale = FALSE, alpha = alpha)
    od_ratio <- if (pc@cutoff.od > 0) pc@od / pc@cutoff.od else 0
    list(score = pmax(pc@sd / pc@cutoff.sd, od_ratio),
         sd = pc@sd, od = pc@od, k = pc@k,
         pcs = pc@scores[, 1:min(2, pc@k), drop = FALSE])
}

#' Returns a tibble with one row per pool in `d` (pools with all metrics).
robust_pca_outliers <- function(d, alpha, pca_metrics = metrics) {
    set.seed(RULE_SEED)
    Z <- robust_z(as.matrix(d[, pca_metrics]))
    n <- nrow(Z)
    obs <- robpca_scores(Z, alpha)
    # Parametric bootstrap under a clean multivariate normal model:
    R <- cov2cor(rrcov::getCov(rrcov::CovMcd(Z)))
    boot <- map(seq_len(N_BOOT), \(i) {
        Zb <- robust_z(MASS::mvrnorm(n, rep(0, ncol(Z)), R))
        robpca_scores(Zb, alpha)$score
    }) |>
        unlist()
    threshold <- quantile(boot, 1 - FP_TARGET, names = FALSE)
    tibble(pool = d$pool, k = obs$k, score = obs$score,
           sd_dist = obs$sd, od_dist = obs$od, threshold = threshold,
           pc1 = obs$pcs[, 1],
           pc2 = if (ncol(obs$pcs) > 1) obs$pcs[, 2] else NA_real_)
}

#' Steps 0-5 on the pools in `d`.
apply_rule <- function(d, max_missing, alpha, pca_metrics = metrics) {
    d <- d |>
        mutate(missing_metric = !complete.cases(pick(all_of(metrics))),
               over_missing = !missing_metric & prop_missing > max_missing)
    pca <- robust_pca_outliers(filter(d, !missing_metric, !over_missing),
                               alpha, pca_metrics)
    d |>
        left_join(pca, by = "pool") |>
        mutate(outlier = !missing_metric & !over_missing & score > threshold,
               include = !missing_metric & !over_missing & !outlier,
               exclude_reasons = case_when(
                   missing_metric ~ "missing metric",
                   over_missing ~ sprintf("prop_missing > %.2f", max_missing),
                   outlier ~ "robust PCA outlier",
                   TRUE ~ ""))
}

time_df <- filter(qc_df, set == "time")
qc_df <- bind_rows(apply_rule(time_df, MAX_MISSING, ALPHA) |>
                       mutate(include_rule = include,
                              override = pool %in% OVERRIDE_INCLUDE,
                              include = include | override,
                              exclude_reasons = if_else(
                                  override,
                                  paste(exclude_reasons,
                                        "(included by owner's override)"),
                                  exclude_reasons)),
                   qc_df |>
                       filter(set == "space") |>
                       mutate(include = TRUE, include_rule = TRUE,
                              exclude_reasons = "no QC (spatial pool)"))


#' Sensitivity (time pools only): the first-run rule, other missing-data
#' limits (2 x and + 0.25 over the median pool's `prop_missing`), and the
#' PCA without `prop_missing` (left to the limit alone).
med_missing <- median(time_df$prop_missing, na.rm = TRUE)
no_missing <- setdiff(metrics, "prop_missing")
sens_rules <- tribble(
    ~rule,                         ~max_missing,         ~alpha, ~pca_metrics,
    "limit 0.5, alpha 0.5 (used)", MAX_MISSING,          ALPHA,  metrics,
    "first run: no limit, alpha 0.75", Inf,              0.75,   metrics,
    "limit 0.5, alpha 0.75",       MAX_MISSING,          0.75,   metrics,
    "limit 2 x median, alpha 0.5", 2 * med_missing,      ALPHA,  metrics,
    "limit median + 0.25, alpha 0.5", med_missing + 0.25, ALPHA, metrics,
    "limit 0.5, alpha 0.5, PCA without prop_missing", MAX_MISSING, ALPHA,
    no_missing)
sens_df <- sens_rules |>
    pmap(\(rule, max_missing, alpha, pca_metrics) {
        apply_rule(time_df, max_missing, alpha, pca_metrics) |>
            transmute(pool, prop_missing, rule = rule, include)
    }) |>
    list_rbind() |>
    pivot_wider(names_from = rule, values_from = include)
write_csv(sens_df, "_data/pool-qc-sensitivity.csv")
cat(sprintf("Median time-pool prop_missing: %.3f\n", med_missing))
sens_df |>
    filter(apply(select(sens_df, -pool, -prop_missing), 1,
                 \(r) n_distinct(r) > 1)) |>
    arrange(desc(prop_missing)) |>
    print(width = Inf)
sens_df |>
    summarize(across(-c(pool, prop_missing), sum)) |>
    pivot_longer(everything(), names_to = "rule", values_to = "n_included") |>
    print()



# =============================================================================*
# Plots ----
# =============================================================================*

#' Metrics by pool, robust PC scores, and the outlier map (SD vs. OD with
#' cutoffs), colored by inclusion.
plot_set <- function(s) {
    d <- qc_df |> filter(set == s)
    p1 <- d |>
        select(pool, include, all_of(metrics)) |>
        pivot_longer(all_of(metrics), names_to = "metric") |>
        ggplot(aes(value, reorder(pool, value), color = include)) +
        geom_point() +
        facet_wrap(~ metric, scales = "free_x", nrow = 1) +
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "firebrick")) +
        theme_classic() +
        theme(axis.title = element_blank(),
              axis.text.y = element_text(size = 6))
    p2 <- d |>
        filter(!is.na(score)) |>
        ggplot(aes(pc1, pc2, label = pool, color = include)) +
        geom_text(size = 2) +
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "firebrick"),
                           guide = "none") +
        theme_classic() +
        labs(x = "Robust PC1", y = "Robust PC2")
    p3 <- d |>
        filter(!is.na(score)) |>
        ggplot(aes(score, reorder(pool, score), color = include)) +
        geom_vline(aes(xintercept = threshold), linetype = 2,
                   color = "gray50") +
        geom_point() +
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "firebrick"),
                           guide = "none") +
        theme_classic() +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_text(size = 6)) +
        labs(x = "Outlier score (dashed: calibrated threshold)",
             caption = sprintf("ROBPCA with k = %i PCs", d$k[1]))
    p1 / (p2 | p3) + plot_annotation(title = paste("Pool QC:", s))
}

ggsave("_plots/pool-qc-time.pdf", plot_set("time"), width = 10, height = 9)



# =============================================================================*
# Outputs ----
# =============================================================================*

write_csv(qc_df, "_data/pool-qc.csv")
print(select(qc_df, set, pool, all_of(metrics), score, threshold,
             include_rule, include,
             exclude_reasons), n = Inf, width = Inf)

for (s in c("time", "space")) {
    qc_df |>
        filter(set == s, include) |>
        getElement("pool") |>
        write_lines(sprintf("_bash/samples/%s-pools.txt", s))
}
qc_df |>
    group_by(set) |>
    summarize(n_candidates = n(), n_included_rule = sum(include_rule),
              n_included = sum(include),
              k = first(na.omit(k))) |>
    print()
