#'
#' Why does SN-1-95_S59 diverge less from the other pools than its BayPass
#' Ω implies? (Owner, 2026-11-19.) `_data/baypass-pod-diagnose.txt`: its
#' mean pairwise divergence from reads is 0.32x the divergence its Ω gives,
#' vs. 0.80-1.11x for every other pool; simulated data at its depth give
#' 1.02x, so the method works at that depth. SN-1-95 was flagged by the
#' pool QC rule and kept by owner override (2026-11-02); median depth 18.
#'
#' Data: mp050 BayPass inputs (3 subsamples; reference read counts and
#' depth per pool; 0 = masked), Ω of mp050 subsample 1. Divergence:
#' moment estimate of Δ_jl = Ω_jj + Ω_ll - 2 Ω_jl from reads, with
#' sampling variance removed (as in `_R/baypass-pod-diagnose.R`).
#'
#'   1. Δ from reads / Δ from Ω for each pair with SN-1-95 (is it low with
#'      every pool, or with some?).
#'   2. The same by SN-1-95's depth at the SNP (if the sampling correction
#'      were wrong at low depth, the ratio would change with depth).
#'   3. A deep pool (SN-1-97_S61, median depth ~100) thinned at random to
#'      SN-1-95's depth at each SNP: does low depth alone lower the ratio?
#'   4. Mixture fit: each pool's frequencies regressed on the other 17
#'      pools' (non-negative least squares, SNPs with data in all pools);
#'      residual variance beyond sampling noise, relative to the
#'      frequency variance. A pool that is partly a mix of other pools'
#'      DNA is unusually well predicted, and its weights show which pools.
#'
#' Writes `_data/sn-1-95-check.txt`.
#'

library(tidyverse)

#' Analyses the first 2026-11 rerun (18 pools; archived):
ARCHIVE_18POOLS <- TRUE
source("_R/data-dir.R")

out_file <- "_data/sn-1-95-check.txt"
sink(out_file, split = TRUE)

set.seed(1306524479)
pools <- read_lines("_bash/samples/time-pools-18.txt")
J <- length(pools)
focal <- "SN-1-95_S59"
deep <- "SN-1-97_S61"
jf <- which(pools == focal)
jd <- which(pools == deep)



# =============================================================================*
# ---- data ----
# =============================================================================*


tmp <- tempfile()
untar(tany_dir("baypass-inputs", "baypass-inputs_mp050_3.tar.gz"), exdir = tmp)
d <- file.path(tmp, "baypass-inputs_mp050")
g <- map(1:3, \(k) data.table::fread(file.path(d, paste0("tany.genobaypass.sub", k)),
                                     header = FALSE) |> as.matrix()) |>
    do.call(what = rbind)
hap <- read_lines(file.path(d, "tany.poolsize")) |>
    str_split_1(" ") |>
    as.integer()
unlink(tmp, recursive = TRUE)
ref <- g[, seq(1, ncol(g), 2)]
dep <- ref + g[, seq(2, ncol(g), 2)]
rm(g)

untar(tany_dir("baypass-init", "baypass-init-mp050-sub1.tar.gz"),
      files = "baypass-init-mp050-sub1/tany-init-mp050-sub1_mat_omega.out",
      exdir = tmp)
om <- as.matrix(read.table(file.path(tmp, "baypass-init-mp050-sub1",
                                     "tany-init-mp050-sub1_mat_omega.out")))
d_om <- outer(diag(om), diag(om), "+") - 2 * om
cat("SNPs:", nrow(ref), "; pools:", J, "\n")



# =============================================================================*
# ---- 1-3. pairwise divergence ----
# =============================================================================*


#' Δ between pools a and b over SNPs `rows`, as in
#' `_R/baypass-pod-diagnose.R`, with e = mean frequency over pools with
#' data.
p_all <- ifelse(dep > 0, ref / pmax(dep, 1), NA_real_)
e_all <- rowMeans(p_all, na.rm = TRUE)
v_all <- e_all * (1 - e_all)
cc <- function(dp, h) 1 / pmax(dp, 1) + 1 / h - 1 / (pmax(dp, 1) * h)
delta <- function(pa, da, ha, pb, db, hb, rows) {
    b <- rows & da > 0 & db > 0
    sum((pa[b] - pb[b])^2 - v_all[b] * (cc(da[b], ha) + cc(db[b], hb))) /
        sum(v_all[b])
}
all_rows <- rep(TRUE, nrow(ref))

pair_tab <- map(setdiff(seq_len(J), jf), \(j) {
    tibble(partner = pools[j],
           ratio = delta(p_all[, jf], dep[, jf], hap[jf], p_all[, j], dep[, j],
                         hap[j], all_rows) / d_om[jf, j])
}) |>
    list_rbind()
cat("\n== 1. Delta(reads) / Delta(Omega), SN-1-95 with each pool\n")
print(as.data.frame(pair_tab |> mutate(ratio = round(ratio, 3)) |> arrange(ratio)))
cat("Mean over partners:", round(mean(pair_tab$ratio), 3), "\n")
#' For comparison: each other pool's mean ratio over its partners.
others <- map_dbl(setdiff(seq_len(J), jf), \(i) {
    mean(map_dbl(setdiff(seq_len(J), i), \(j) {
        delta(p_all[, i], dep[, i], hap[i], p_all[, j], dep[, j], hap[j],
              all_rows) / d_om[i, j]
    }))
})
cat("Other pools' mean ratios: range", paste(round(range(others), 3),
                                             collapse = "-"), "\n")

cat("\n== 2. By SN-1-95 depth at the SNP (mean ratio over partners)\n")
bins <- cut(dep[, jf], c(0, 12, 18, 25, Inf), right = TRUE)
map(levels(bins), \(lv) {
    rows <- !is.na(bins) & bins == lv
    tibble(depth = lv, snps = sum(rows),
           ratio = mean(map_dbl(setdiff(seq_len(J), jf), \(j) {
               delta(p_all[, jf], dep[, jf], hap[jf], p_all[, j], dep[, j],
                     hap[j], rows) / d_om[jf, j]
           })))
}) |>
    list_rbind() |>
    mutate(ratio = round(ratio, 3)) |>
    as.data.frame() |>
    print()

cat("\n== 3.", deep, "thinned to SN-1-95's depth at each SNP\n")
both <- dep[, jd] > 0 & dep[, jf] > 0
td <- pmin(dep[, jf], dep[, jd])
tr <- ifelse(both, rhyper(nrow(ref), ref[, jd], dep[, jd] - ref[, jd], td), 0L)
p_thin <- ifelse(both & td > 0, tr / pmax(td, 1), NA_real_)
ratio_thin <- map_dbl(setdiff(seq_len(J), c(jd, jf)), \(j) {
    delta(p_thin, ifelse(both, td, 0L), hap[jd], p_all[, j], dep[, j], hap[j],
          all_rows) / d_om[jd, j]
})
ratio_full <- map_dbl(setdiff(seq_len(J), c(jd, jf)), \(j) {
    delta(p_all[, jd], dep[, jd], hap[jd], p_all[, j], dep[, j], hap[j],
          both) / d_om[jd, j]
})
cat("Mean ratio over partners: full depth", round(mean(ratio_full), 3),
    "; thinned", round(mean(ratio_thin), 3), "\n")



# =============================================================================*
# ---- 4. mixture fit ----
# =============================================================================*


#' SNPs with data in all pools and depth >= 10 everywhere.
ok <- rowSums(dep >= 10) == J
P <- p_all[ok, ]
D <- dep[ok, ]
cat("\n== 4. Mixture fit (", sum(ok), " SNPs with depth >= 10 in all pools)\n",
    sep = "")
mix <- map(seq_len(J), \(t) {
    X <- P[, -t]
    f <- nnls::nnls(X, P[, t])
    w <- setNames(f$x, pools[-t])
    fit <- as.vector(X %*% f$x)
    samp <- mean(v_all[ok] * cc(D[, t], hap[t]))
    tibble(pool = pools[t], sum_w = sum(w),
           excess_resid = (mean((P[, t] - fit)^2) - samp) / mean(v_all[ok]),
           top_weights = paste(names(sort(w, decreasing = TRUE))[1:3],
                               round(sort(w, decreasing = TRUE)[1:3], 2),
                               collapse = "; "))
}) |>
    list_rbind()
cat("(excess_resid: residual variance beyond sampling noise, over mean",
    "e(1 - e); low = well predicted by the other pools)\n")
print(as.data.frame(mix |> mutate(across(c(sum_w, excess_resid),
                                         \(x) round(x, 4))) |>
                        arrange(excess_resid)))

#' Raw divergences behind the ratios for SN-1-95's closest partners, and
#' each partner's typical divergence from the other pools (Δ from reads).
cat("\n== 5. Delta from reads and from Omega, SN-1-95 vs. its closest",
    "partners\n")
d_reads_95 <- map_dbl(seq_len(J), \(j) {
    if (j == jf) return(NA_real_)
    delta(p_all[, jf], dep[, jf], hap[jf], p_all[, j], dep[, j], hap[j], all_rows)
})
typ <- map_dbl(seq_len(J), \(j) {
    median(map_dbl(setdiff(seq_len(J), c(j, jf)), \(k) {
        delta(p_all[, j], dep[, j], hap[j], p_all[, k], dep[, k], hap[k], all_rows)
    }))
})
tibble(partner = pools, delta_reads = d_reads_95, delta_omega = d_om[jf, ],
       partner_typical_delta = typ) |>
    filter(partner != focal) |>
    arrange(delta_reads) |>
    mutate(across(where(is.double), \(x) signif(x, 3))) |>
    as.data.frame() |>
    print()
cat("(partner_typical_delta: median Delta from reads between that pool and",
    "the other pools, SN-1-95 excluded)\n")

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
