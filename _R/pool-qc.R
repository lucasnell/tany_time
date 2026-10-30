
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
#' Rule (decided 2026-10-10, before any QC metrics existed): within each
#' sample set (time, space), a pool is excluded if it is an outlier in a
#' robust PCA of the QC metrics, or if any metric is missing.
#'
#'   1. Standardize each metric robustly: (x - median) / MAD.
#'   2. ROBPCA (`rrcov::PcaHubert`; Hubert et al. 2005,
#'      doi:10.1198/004017004000000563) with its default number of PCs and
#'      default (0.975) cutoffs for score distance (SD; within the PCs) and
#'      orthogonal distance (OD; off the PCs).
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
#' Steps:
#'   1. Download and extract `pool-qc.tar` into `qc_dir` (below).
#'   2. Run this script. It writes `_bash/samples/{time,space}-pools.txt`,
#'      used by `02-allele-freqs/04-combine-syncs`.
#'
#' Outputs: `_data/pool-qc.csv`, `_plots/pool-qc-{time,space}.pdf`, and the
#' pool lists above.
#'

library(tidyverse)
library(patchwork)

qc_dir <- "~/_data/_tany-time/pool-qc"

#' QC metrics used in the rule (definitions in `03-pool-qc.sh` and
#' `pool-qc-files.py`). The mapping rate is `target_rate`: reads assigned to
#' T. gracilentus contigs (MAPQ >= 20, primary) under competitive mapping,
#' over all reads, excluding reads on a mitochondrial-like stretch of
#' CONTIG37 (`mitolike_rate`; decided 2026-10-30, before any QC metrics
#' existed, because they are mitochondrial reads).
metrics <- c("prop_missing", "n_private", "picard_dup", "fastp_gc",
             "target_rate")

#' Rule parameters (fixed before seeing any metrics; see above):
FP_TARGET <- 0.025     # per-pool false-positive rate under the clean model
N_BOOT <- 2000         # simulated data sets for calibration
RULE_SEED <- 1702635893 # MCD and bootstrap use random numbers



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
robpca_scores <- function(Z) {
    pc <- rrcov::PcaHubert(Z, k = 0, scale = FALSE)
    od_ratio <- if (pc@cutoff.od > 0) pc@od / pc@cutoff.od else 0
    list(score = pmax(pc@sd / pc@cutoff.sd, od_ratio),
         sd = pc@sd, od = pc@od, k = pc@k,
         pcs = pc@scores[, 1:min(2, pc@k), drop = FALSE])
}

#' Returns a tibble with one row per pool in `d`. Pools with any missing
#' metric get NA (excluded later).
robust_pca_outliers <- function(d) {
    set.seed(RULE_SEED)
    complete <- complete.cases(d[, metrics])
    Z <- robust_z(as.matrix(d[complete, metrics]))
    n <- nrow(Z)
    obs <- robpca_scores(Z)
    # Parametric bootstrap under a clean multivariate normal model:
    R <- cov2cor(rrcov::getCov(rrcov::CovMcd(Z)))
    boot <- map(seq_len(N_BOOT), \(i) {
        Zb <- robust_z(MASS::mvrnorm(n, rep(0, ncol(Z)), R))
        robpca_scores(Zb)$score
    }) |>
        unlist()
    threshold <- quantile(boot, 1 - FP_TARGET, names = FALSE)
    out <- tibble(pool = d$pool, k = obs$k, score = NA_real_,
                  sd_dist = NA_real_, od_dist = NA_real_,
                  threshold = threshold, pc1 = NA_real_, pc2 = NA_real_)
    out$score[complete] <- obs$score
    out$sd_dist[complete] <- obs$sd
    out$od_dist[complete] <- obs$od
    out$pc1[complete] <- obs$pcs[, 1]
    if (ncol(obs$pcs) > 1) out$pc2[complete] <- obs$pcs[, 2]
    return(out)
}

qc_df <- qc_df |>
    group_by(set) |>
    group_modify(\(d, g) left_join(d, robust_pca_outliers(d), by = "pool")) |>
    ungroup() |>
    mutate(missing_metric = !complete.cases(pick(all_of(metrics))),
           outlier = !missing_metric & score > threshold,
           include = !missing_metric & !outlier,
           exclude_reasons = case_when(missing_metric ~ "missing metric",
                                       outlier ~ "robust PCA outlier",
                                       TRUE ~ ""))



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
        filter(!missing_metric) |>
        ggplot(aes(pc1, pc2, label = pool, color = include)) +
        geom_text(size = 2) +
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "firebrick"),
                           guide = "none") +
        theme_classic() +
        labs(x = "Robust PC1", y = "Robust PC2")
    p3 <- d |>
        filter(!missing_metric) |>
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

for (s in c("time", "space")) {
    ggsave(sprintf("_plots/pool-qc-%s.pdf", s), plot_set(s),
           width = 10, height = 9)
}



# =============================================================================*
# Outputs ----
# =============================================================================*

write_csv(qc_df, "_data/pool-qc.csv")
print(select(qc_df, set, pool, all_of(metrics), score, threshold, include,
             exclude_reasons), n = Inf, width = Inf)

for (s in c("time", "space")) {
    qc_df |>
        filter(set == s, include) |>
        getElement("pool") |>
        write_lines(sprintf("_bash/samples/%s-pools.txt", s))
}
qc_df |>
    group_by(set) |>
    summarize(n_candidates = n(), n_included = sum(include), k = k[1]) |>
    print()
