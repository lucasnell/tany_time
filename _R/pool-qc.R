
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
#' Rule (decided 2026-10-10, before any QC metrics existed): within each
#' sample set (time, space), a pool is excluded if it is an outlier in a
#' robust PCA of the QC metrics, or if any metric is missing.
#'
#'   1. Standardize each metric robustly: (x - median) / MAD.
#'   2. Robust covariance by minimum covariance determinant (MCD;
#'      `MASS::cov.rob`), so outlying pools don't shape the PCs.
#'   3. PCs = eigenvectors of that covariance; keep the first k explaining
#'      >= 80% of variance (k <= number of metrics - 1).
#'   4. Score distance (SD): Mahalanobis distance within the k PCs; cutoff
#'      sqrt(chi-squared_k quantile 0.975).
#'      Orthogonal distance (OD): distance from the k-PC subspace; cutoff
#'      (median(OD^(2/3)) + MAD(OD^(2/3)) * z_0.975)^(3/2).
#'      These are the ROBPCA cutoffs of Hubert et al. 2005
#'      (doi:10.1198/004017004000000563).
#'   5. Outlier if SD or OD exceeds its cutoff.
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

#' QC metrics used in the rule (see `03-pool-qc.sh` for definitions):
metrics <- c("prop_missing", "n_private", "picard_dup", "fastp_gc",
             "map_rate")

#' Rule parameters (fixed before seeing any metrics; see above):
PCA_PROP_VAR <- 0.8    # proportion of variance the kept PCs must explain
CUTOFF_Q <- 0.975      # quantile for SD and OD cutoffs
MCD_SEED <- 1702635893 # MCD uses random subsets; fixed for reproducibility

#' NOT YET CALIBRATED. `_R/pool-qc-rule-sim.R` shows this implementation
#' (`MASS::cov.rob` without small-sample corrections) excludes 2–7 clean
#' pools per set on simulated clean data. The script stops before writing
#' pool lists until the rule is fixed (see TODO.md).
RULE_CALIBRATED <- FALSE



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




# =============================================================================*
# Robust PCA outlier rule ----
# =============================================================================*

#' Returns a tibble with one row per pool in `d`: SD, OD, cutoffs, k, and
#' PC1/PC2 scores. Pools with any missing metric get NA (excluded later).
robust_pca_outliers <- function(d) {
    complete <- complete.cases(d[, metrics])
    X <- as.matrix(d[complete, metrics])
    Z <- apply(X, 2, \(x) {
        s <- mad(x)
        if (s == 0) s <- sd(x)
        if (s == 0) stop("metric with no variation")
        (x - median(x)) / s
    })
    set.seed(MCD_SEED)
    mcd <- MASS::cov.rob(Z, method = "mcd")
    eig <- eigen(mcd$cov, symmetric = TRUE)
    p <- ncol(Z)
    k <- min(which(cumsum(eig$values) / sum(eig$values) >= PCA_PROP_VAR),
             p - 1L)
    Zc <- sweep(Z, 2, mcd$center)
    V <- eig$vectors[, 1:k, drop = FALSE]
    scores <- Zc %*% V
    sd_dist <- sqrt(rowSums(sweep(scores^2, 2, eig$values[1:k], "/")))
    od_dist <- sqrt(rowSums((Zc - scores %*% t(V))^2))
    sd_cut <- sqrt(qchisq(CUTOFF_Q, k))
    od23 <- od_dist^(2/3)
    od_cut <- (median(od23) + mad(od23) * qnorm(CUTOFF_Q))^(3/2)
    out <- tibble(pool = d$pool, k = k, sd_dist = NA_real_,
                  od_dist = NA_real_, sd_cut = sd_cut, od_cut = od_cut,
                  pc1 = NA_real_, pc2 = NA_real_,
                  prop_var_k = sum(eig$values[1:k]) / sum(eig$values))
    out$sd_dist[complete] <- sd_dist
    out$od_dist[complete] <- od_dist
    out$pc1[complete] <- scores[, 1]
    if (k > 1) out$pc2[complete] <- scores[, 2]
    return(out)
}

qc_df <- qc_df |>
    group_by(set) |>
    group_modify(\(d, g) left_join(d, robust_pca_outliers(d), by = "pool")) |>
    ungroup() |>
    mutate(missing_metric = !complete.cases(pick(all_of(metrics))),
           outlier_sd = !missing_metric & sd_dist > sd_cut,
           outlier_od = !missing_metric & od_dist > od_cut,
           include = !missing_metric & !outlier_sd & !outlier_od,
           exclude_reasons = pmap_chr(
               list(missing_metric, outlier_sd, outlier_od),
               \(m, s, o) paste(c("missing metric", "score distance",
                                  "orthogonal distance")[c(m, s, o)],
                                collapse = ";")))



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
        ggplot(aes(sd_dist, od_dist, label = pool, color = include)) +
        geom_vline(aes(xintercept = sd_cut), linetype = 2, color = "gray50") +
        geom_hline(aes(yintercept = od_cut), linetype = 2, color = "gray50") +
        geom_text(size = 2) +
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "firebrick"),
                           guide = "none") +
        theme_classic() +
        labs(x = "Score distance", y = "Orthogonal distance",
             caption = sprintf("k = %i PCs (%.0f%% of robust variance)",
                               d$k[1], 100 * d$prop_var_k[1]))
    p1 / (p2 | p3) + plot_annotation(title = paste("Pool QC:", s))
}

for (s in c("time", "space")) {
    ggsave(sprintf("_plots/pool-qc-%s.pdf", s), plot_set(s),
           width = 10, height = 9)
}



# =============================================================================*
# Outputs ----
# =============================================================================*

if (!RULE_CALIBRATED) {
    stop("Outlier rule not calibrated; no outputs written (see TODO.md)")
}

write_csv(qc_df, "_data/pool-qc.csv")
print(select(qc_df, set, pool, all_of(metrics), sd_dist, od_dist, include,
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
