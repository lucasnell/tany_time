#'
#' Why does Omega differ among missing-data settings? (2026-11-05)
#'
#' `baypass-check-omegas.R` found FMD between mp050 and mp100 Omegas
#' (0.56-0.64) larger than between subsamples within a setting
#' (0.22-0.33). The settings differ only in the minimum proportion of pools
#' with data (after masking) for a SNP to be kept (0.5, 0.75, 1), so
#' mp050 = mp100-type SNPs (data in all pools) + SNPs with missing pools.
#' Possible causes:
#'   1. scale: Omega multiplied by a constant (FMD is sensitive to scale);
#'   2. particular pools: rows/columns of a few pools differ;
#'   3. the SNPs with missing pools differ in their covariance (which SNPs
#'      are included), or missing data itself biases the estimate.
#'
#' For 3, BayPass can't be rerun locally (10-30 h per run), so this uses a
#' moment estimate of pairwise divergence (defined below), computed on (a) all mp050
#' SNPs, (b) mp050 SNPs with data in all pools, (c) mp050 SNPs with missing
#' pools, (d) mp100 SNPs, and (e) the complete SNPs of (b) with the
#' missing-data pattern of (c) imposed on them.
#'
#' Reads `tany_dir("baypass-init")` and `tany_dir("baypass-inputs")`.
#' Writes `_data/baypass-omega-diffs.txt` (all printed output).
#'

library(tidyverse)
library(geigen) # geigen function inside fmd.dist

source("_R/data-dir.R")

out_file <- "_data/baypass-omega-diffs.txt"
sink(out_file, split = TRUE)

mp_tags <- c("mp050", "mp075", "mp100")
#' Pool order of the BayPass inputs (see `_R/baypass-covariates.R`):
pools <- read_lines("_bash/samples/time-pools.txt")
n_pools <- length(pools)

fmd.dist <- function(mat1, mat2) {
    return(sqrt(sum(log(geigen(mat1,mat2)$values)**2)))
}



# =============================================================================*
# ---- read Omegas and inputs ----
# =============================================================================*


read_omega <- function(tag, k) {
    nm <- sprintf("%s-sub%i", tag, k)
    fn <- sprintf("baypass-init-%s/tany-init-%s_mat_omega.out", nm, nm)
    untar(tany_dir("baypass-init", sprintf("baypass-init-%s.tar.gz", nm)),
          files = fn, exdir = tempdir())
    file.path(tempdir(), fn) |>
        read_table(col_names = FALSE, col_types = cols()) |>
        as.matrix() |>
        unname()
}

#' Reference and alternative read counts (SNPs x pools) and haploid pool
#' sizes for one setting and subsample.
read_counts <- function(tag, k) {
    d <- sprintf("baypass-inputs_%s", tag)
    fns <- file.path(d, c(paste0("tany.genobaypass.sub", k), "tany.poolsize"))
    untar(tany_dir("baypass-inputs", sprintf("baypass-inputs_%s_3.tar.gz", tag)),
          files = fns, exdir = tempdir())
    g <- file.path(tempdir(), fns[1]) |>
        data.table::fread(header = FALSE) |>
        as.matrix() |>
        unname()
    h <- file.path(tempdir(), fns[2]) |>
        read_lines() |>
        str_split_1(" ") |>
        as.integer()
    list(ref = g[, seq(1, ncol(g), 2)], alt = g[, seq(2, ncol(g), 2)],
         hap = h)
}

omegas <- crossing(tag = mp_tags, k = 1:3) |>
    mutate(om = map2(tag, k, read_omega))

counts <- crossing(tag = mp_tags, k = 1:3) |>
    mutate(cnt = map2(tag, k, read_counts))



# =============================================================================*
# ---- 1. scale ----
# =============================================================================*


#' Generalized eigenvalues of (Omega_a, Omega_b), matched subsamples (sub k
#' vs. sub k). FMD = sqrt(sum(log(lambda)^2)). If Omega_a = c * Omega_b,
#' all lambda = c. `scaled_fmd` is FMD after dividing each Omega by its
#' mean diagonal.
om_of <- function(tag, k) omegas$om[[which(omegas$tag == tag & omegas$k == k)]]

scale_df <- crossing(a = mp_tags, b = mp_tags, k = 1:3) |>
    filter(a < b) |>
    mutate(lam = pmap(list(a, b, k), \(a, b, k) {
        sort(geigen(om_of(a, k), om_of(b, k))$values)
    }),
    fmd = map_dbl(lam, \(x) sqrt(sum(log(x)^2))),
    mean_log_lam = map_dbl(lam, \(x) mean(log(x))),
    lam_min = map_dbl(lam, min),
    lam_max = map_dbl(lam, max),
    trace_ratio = pmap_dbl(list(a, b, k), \(a, b, k) {
        mean(diag(om_of(a, k))) / mean(diag(om_of(b, k)))
    }),
    scaled_fmd = pmap_dbl(list(a, b, k), \(a, b, k) {
        oa <- om_of(a, k)
        ob <- om_of(b, k)
        fmd.dist(oa / mean(diag(oa)), ob / mean(diag(ob)))
    }))

cat("\n== 1. Scale: generalized eigenvalues, matched subsamples\n")
scale_df |>
    select(-lam) |>
    print(n = Inf, width = Inf)

cat("\nEigenvalues (sorted), sub1 of each pair:\n")
scale_df |>
    filter(k == 1) |>
    mutate(lam = map_chr(lam, \(x) paste(round(x, 2), collapse = " "))) |>
    select(a, b, lam) |>
    as.data.frame() |>
    print()

#' Pools loading on the most extreme generalized eigenvectors (mp050 vs.
#' mp100): squared loadings of the eigenvector scaled to sum to 1, mean
#' over subsamples. These are the directions FMD mostly comes from.
ext_df <- map(1:3, \(k) {
    ge <- geigen(om_of("mp050", k), om_of("mp100", k))
    v_min <- ge$vectors[, which.min(ge$values)]
    v_max <- ge$vectors[, which.max(ge$values)]
    tibble(pool = pools, k = k,
           load_min = v_min^2 / sum(v_min^2),
           load_max = v_max^2 / sum(v_max^2))
}) |>
    list_rbind() |>
    group_by(pool) |>
    summarize(across(c(load_min, load_max), mean)) |>
    arrange(desc(load_min))
cat("\nSquared loadings on the smallest / largest eigenvector (mp050 vs. mp100):\n")
ext_df |>
    mutate(across(where(is.double), \(x) round(x, 3))) |>
    print(n = Inf)



# =============================================================================*
# ---- 2. pools ----
# =============================================================================*


#' Per pool: mean Omega diagonal per setting (over subsamples), proportion
#' of SNPs with no reads (missing), and median read depth where present.
pool_df <- counts |>
    mutate(df = map(cnt, \(x) {
        dep <- x$ref + x$alt
        tibble(pool = pools, hap = x$hap,
               missing = colMeans(dep == 0),
               depth = apply(dep, 2, \(d) median(d[d > 0])))
    })) |>
    select(tag, k, df) |>
    unnest(df) |>
    left_join(omegas |>
                  mutate(df = map(om, \(x) tibble(pool = pools,
                                                  omega_diag = diag(x)))) |>
                  select(tag, k, df) |>
                  unnest(df),
              by = c("tag", "k", "pool")) |>
    group_by(tag, pool, hap) |>
    summarize(across(c(missing, depth, omega_diag), mean), .groups = "drop") |>
    pivot_wider(names_from = tag, values_from = c(missing, depth, omega_diag))

#' FMD between mp050 and mp100 (matched subsamples, mean over k) after
#' dropping one pool from both matrices. A pool whose removal lowers FMD
#' most contributes most to the difference.
fmd_all <- mean(map_dbl(1:3, \(k) fmd.dist(om_of("mp050", k), om_of("mp100", k))))
drop_df <- tibble(pool = pools, i = seq_along(pools)) |>
    mutate(fmd_drop = map_dbl(i, \(i) {
        mean(map_dbl(1:3, \(k) {
            fmd.dist(om_of("mp050", k)[-i, -i], om_of("mp100", k)[-i, -i])
        }))
    })) |>
    select(-i)

cat("\n== 2. Pools (means over 3 subsamples)\n")
cat(sprintf("FMD mp050 vs. mp100 with all pools: %.3f\n", fmd_all))
pool_df |>
    left_join(drop_df, by = "pool") |>
    mutate(diag_ratio = omega_diag_mp050 / omega_diag_mp100) |>
    select(pool, hap, missing_mp050, missing_mp075, depth_mp050,
           depth_mp100, omega_diag_mp050, omega_diag_mp100, diag_ratio,
           fmd_drop) |>
    arrange(fmd_drop) |>
    mutate(across(where(is.double), \(x) signif(x, 3))) |>
    print(n = Inf, width = Inf)



# =============================================================================*
# ---- 3. which SNPs ----
# =============================================================================*


#' Pairwise divergence. Omega itself can't be estimated by moments
#' without centering each SNP on a mean over pools, and that mean depends
#' on which pools have data, so missing data alone would change it. The
#' pairwise divergence Delta_jl = Omega_jj + Omega_ll - 2 Omega_jl doesn't
#' need a center, and uses only pools j and l. For SNP i, p_ij = ref reads
#' / depth, e_i = mean of p_ij over pools with data, and the sampling
#' variance of p_ij (reads and pool) is about e_i (1 - e_i) c_ij, with
#' c_ij = 1/d_ij + 1/h_j - 1/(d_ij h_j) (d = depth, h = haploid pool
#' size). Over SNPs where both pools have data:
#'   Delta_jl = sum_i [(p_ij - p_il)^2 - e_i (1 - e_i) (c_ij + c_il)]
#'              / sum_i e_i (1 - e_i).
#' This is not BayPass's estimator; it compares SNP sets with one method.
#' The same Delta computed from each BayPass Omega checks that the method
#' sees the mp050 vs. mp100 difference.
#' Optional SNP weights `w` (section 5) multiply each SNP's terms in both
#' sums.
pair_delta <- function(ref, alt, hap, w = rep(1, nrow(ref))) {
    dep <- ref + alt
    has <- dep > 0
    p <- ifelse(has, ref / pmax(dep, 1), NA_real_)
    e <- rowMeans(p, na.rm = TRUE)
    v <- e * (1 - e)
    cc <- sweep(1 / pmax(dep, 1), 2, 1 / hap, "+") -
        sweep(1 / pmax(dep, 1), 2, hap, "/")
    J <- ncol(dep)
    out <- matrix(0, J, J)
    for (j in 1:(J - 1)) {
        for (l in (j + 1):J) {
            b <- has[, j] & has[, l]
            out[j, l] <- out[l, j] <-
                sum(w[b] * ((p[b, j] - p[b, l])^2 -
                                v[b] * (cc[b, j] + cc[b, l]))) /
                sum(w[b] * v[b])
        }
    }
    out
}
delta_of_omega <- function(om) {
    outer(diag(om), diag(om), "+") - 2 * om
}

sets_one <- function(k) {
    c050 <- counts$cnt[[which(counts$tag == "mp050" & counts$k == k)]]
    c100 <- counts$cnt[[which(counts$tag == "mp100" & counts$k == k)]]
    complete <- rowSums(c050$ref + c050$alt == 0) == 0
    #' (e): complete SNPs with the missing pattern of incomplete SNPs,
    #' taken in order and recycled.
    miss <- (c050$ref + c050$alt == 0)[!complete, , drop = FALSE]
    rc <- c050$ref[complete, ]
    ac <- c050$alt[complete, ]
    idx <- rep_len(seq_len(nrow(miss)), nrow(rc))
    rc_e <- ifelse(miss[idx, ], 0L, rc)
    ac_e <- ifelse(miss[idx, ], 0L, ac)
    list(
        n = c(all = nrow(c050$ref), complete = sum(complete),
              incomplete = sum(!complete), mp100 = nrow(c100$ref)),
        n_missing = table(rowSums(c050$ref + c050$alt == 0)),
        delta = list(
            mp050_all = pair_delta(c050$ref, c050$alt, c050$hap),
            mp050_complete = pair_delta(rc, ac, c050$hap),
            mp050_incomplete = pair_delta(c050$ref[!complete, ],
                                          c050$alt[!complete, ], c050$hap),
            mp100 = pair_delta(c100$ref, c100$alt, c100$hap),
            complete_with_missing = pair_delta(rc_e, ac_e, c050$hap)))
}

mom <- map(1:3, sets_one)

#' Sets are compared by D(A, B) = sqrt(sum((A_jl - B_jl)^2)) /
#' sqrt(sum(B_jl^2)) over pairs j < l. Noise baseline: D between
#' subsamples of the same set.
off <- function(x) x[upper.tri(x)]
dist_off <- function(a, b) sqrt(sum((off(a) - off(b))^2)) / sqrt(sum(off(b)^2))

cat("\n== 3. Pairwise divergence (Delta) by SNP set\n")
cat("SNPs per set (sub1, sub2, sub3):\n")
print(sapply(mom, \(x) x$n))
cat("\nmp050 sub1: SNPs by number of pools missing:\n")
print(mom[[1]]$n_missing)

cat("\nBayPass Omegas converted to Delta, D:\n")
bp_d <- bind_rows(
    crossing(a = mp_tags, b = mp_tags, k = 1:3) |>
        filter(a < b) |>
        mutate(d = pmap_dbl(list(a, b, k), \(a, b, k) {
            dist_off(delta_of_omega(om_of(a, k)), delta_of_omega(om_of(b, k)))
        }),
        pair = paste(a, "vs.", b, "(matched sub)")),
    crossing(a = mp_tags, k1 = 1:3, k2 = 1:3) |>
        filter(k1 < k2) |>
        mutate(d = pmap_dbl(list(a, k1, k2), \(t, i, j) {
            dist_off(delta_of_omega(om_of(t, i)), delta_of_omega(om_of(t, j)))
        }),
        pair = paste(a, "subsamples")))
bp_d |>
    group_by(pair) |>
    summarize(d_min = min(d), d_max = max(d)) |>
    print()

set_names_ <- names(mom[[1]]$delta)
mom_d <- bind_rows(
    map(1:3, \(k) {
        combn(set_names_, 2) |>
            t() |>
            as_tibble(.name_repair = \(x) c("a", "b")) |>
            mutate(d = map2_dbl(a, b, \(i, j) {
                dist_off(mom[[k]]$delta[[i]], mom[[k]]$delta[[j]])
            }))
    }) |>
        list_rbind(),
    crossing(a = set_names_, k1 = 1:3, k2 = 1:3) |>
        filter(k1 < k2) |>
        mutate(b = a,
               d = pmap_dbl(list(a, k1, k2), \(s, i, j) {
                   dist_off(mom[[i]]$delta[[s]], mom[[j]]$delta[[s]])
               })) |>
        select(a, b, d))
cat("\nDelta by SNP set, D (over 3 subsamples; a = b rows are subsample noise):\n")
mom_d |>
    group_by(a, b) |>
    summarize(d_min = min(d), d_max = max(d), .groups = "drop") |>
    arrange(d_min) |>
    print(n = Inf, width = Inf)

cat("\nDelta from reads vs. from BayPass Omega, correlation over pairs:\n")
map(1:3, \(k) {
    tibble(k = k,
           mp050 = cor(off(mom[[k]]$delta$mp050_all),
                       off(delta_of_omega(om_of("mp050", k)))),
           mp100 = cor(off(mom[[k]]$delta$mp100),
                       off(delta_of_omega(om_of("mp100", k)))))
}) |>
    list_rbind() |>
    print()

#' Direction: is divergence higher among incomplete SNPs? Mean ratio over
#' pairs, and per pool (mean over that pool's pairs), mean of subsamples.
ratio_of <- function(a, b) {
    r <- a / b
    diag(r) <- NA
    r
}
cat("\nMean Delta ratio over pairs (mean of subsamples):\n")
tibble(
    incomplete_vs_complete = mean(sapply(1:3, \(k) {
        mean(off(ratio_of(mom[[k]]$delta$mp050_incomplete,
                          mom[[k]]$delta$mp050_complete)))
    })),
    complete_with_missing_vs_complete = mean(sapply(1:3, \(k) {
        mean(off(ratio_of(mom[[k]]$delta$complete_with_missing,
                          mom[[k]]$delta$mp050_complete)))
    })),
    baypass_050_vs_100 = mean(sapply(1:3, \(k) {
        mean(off(ratio_of(delta_of_omega(om_of("mp050", k)),
                          delta_of_omega(om_of("mp100", k)))))
    }))) |>
    print(width = Inf)

cat("\nPer pool, mean Delta ratio over its pairs:\n")
tibble(pool = pools,
       incompl_vs_compl = rowMeans(sapply(1:3, \(k) {
           rowMeans(ratio_of(mom[[k]]$delta$mp050_incomplete,
                             mom[[k]]$delta$mp050_complete), na.rm = TRUE)
       })),
       baypass_050_vs_100 = rowMeans(sapply(1:3, \(k) {
           rowMeans(ratio_of(delta_of_omega(om_of("mp050", k)),
                             delta_of_omega(om_of("mp100", k))), na.rm = TRUE)
       }))) |>
    arrange(desc(incompl_vs_compl)) |>
    mutate(across(where(is.double), \(x) round(x, 3))) |>
    print(n = Inf)



# =============================================================================*
# ---- 4. what the incomplete SNPs are ----
# =============================================================================*


#' A pool is missing at a site when its depth there is < 10 or above its
#' 95th percentile, or the site is within 5 bp of an indel in that pool
#' (`01-snape.sh` partial mask; repeats are masked in all pools). Depth of
#' the masked pool is unknown, so SNPs are described by the pools that do
#' have data: relative depth = mean over those pools of depth / that pool's
#' median depth (over complete SNPs), and the mean frequency e_i (minor
#' allele frequency = min(e, 1 - e)).
cat("\n== 4. Complete vs. incomplete mp050 SNPs (sub1)\n")
c050 <- counts$cnt[[which(counts$tag == "mp050" & counts$k == 1)]]
dep <- c050$ref + c050$alt
has <- dep > 0
complete <- rowSums(!has) == 0
med <- apply(dep[complete, ], 2, median)
rel_dep <- rowSums(sweep(dep, 2, med, "/") * has) / rowSums(has)
e <- rowSums(ifelse(has, c050$ref / pmax(dep, 1), 0)) / rowSums(has)
snp_df <- tibble(n_missing = rowSums(!has), rel_dep = rel_dep,
                 maf = pmin(e, 1 - e)) |>
    mutate(group = case_when(n_missing == 0 ~ "0",
                             n_missing <= 2 ~ "1-2",
                             n_missing <= 5 ~ "3-5",
                             TRUE ~ "6-9"))
snp_df |>
    group_by(group) |>
    summarize(n = n(),
              rel_dep_q10 = quantile(rel_dep, 0.1),
              rel_dep_median = median(rel_dep),
              rel_dep_q90 = quantile(rel_dep, 0.9),
              maf_median = median(maf)) |>
    mutate(across(where(is.double), \(x) round(x, 3))) |>
    print(width = Inf)

#' Which pools are missing in SNPs with 1 missing pool, relative to each
#' pool's overall missing rate.
cat("\nSNPs with exactly 1 missing pool: share by pool\n")
one <- which(rowSums(!has) == 1)
tibble(pool = pools,
       share_of_1missing = colMeans(!has[one, ])) |>
    arrange(desc(share_of_1missing)) |>
    mutate(share_of_1missing = round(share_of_1missing, 3)) |>
    print(n = Inf)

#' Delta of incomplete SNPs split by relative depth (in pools with data):
#' does the difference from complete SNPs come from high-depth sites
#' (likely duplicated or collapsed regions)?
cat("\nIncomplete SNPs split by relative depth: D from complete SNPs\n")
d_compl <- mom[[1]]$delta$mp050_complete
brks <- c(0, 0.8, 1.2, Inf)
inc <- which(!complete)
tibble(bin = cut(rel_dep[inc], brks)) |>
    mutate(i = inc) |>
    group_by(bin) |>
    summarize(n = n(),
              d = dist_off(pair_delta(c050$ref[i, ], c050$alt[i, ], c050$hap),
                           d_compl),
              mean_ratio = mean(off(ratio_of(
                  pair_delta(c050$ref[i, ], c050$alt[i, ], c050$hap),
                  d_compl))),
              .groups = "drop") |>
    print()
cat("Complete SNPs split the same way:\n")
cmp <- which(complete)
tibble(bin = cut(rel_dep[cmp], brks)) |>
    mutate(i = cmp) |>
    group_by(bin) |>
    summarize(n = n(),
              d = dist_off(pair_delta(c050$ref[i, ], c050$alt[i, ], c050$hap),
                           d_compl),
              mean_ratio = mean(off(ratio_of(
                  pair_delta(c050$ref[i, ], c050$alt[i, ], c050$hap),
                  d_compl))),
              .groups = "drop") |>
    print()



# =============================================================================*
# ---- 5. depth-matched comparison ----
# =============================================================================*


#' Does depth explain the difference between incomplete and complete SNPs?
#' Complete SNPs are reweighted to match the incomplete SNPs' distribution
#' of relative depth (20 quantile bins), or of relative depth x minor
#' allele frequency (20 x 5 bins): weight = share of incomplete SNPs in the
#' bin / share of complete SNPs in the bin. If depth (and MAF) explain the
#' difference, D(incomplete, reweighted complete) falls to subsample noise.
#' The same is done with all mp050 SNPs as the target, compared with mp100
#' (the BayPass comparison). n_eff = (sum w)^2 / sum w^2 (effective number
#' of SNPs after weighting).
snp_info <- function(cnt) {
    dep <- cnt$ref + cnt$alt
    has <- dep > 0
    complete <- rowSums(!has) == 0
    med <- apply(dep[complete, ], 2, median)
    rel_dep <- rowSums(sweep(dep, 2, med, "/") * has) / rowSums(has)
    e <- rowSums(ifelse(has, cnt$ref / pmax(dep, 1), 0)) / rowSums(has)
    tibble(complete = complete, rel_dep = rel_dep, maf = pmin(e, 1 - e))
}
match_weights <- function(bin, from, to) {
    #' weights for SNPs `from` so their bin shares equal those of `to`
    p_to <- table(factor(bin[to], levels = levels(bin))) / sum(to)
    p_from <- table(factor(bin[from], levels = levels(bin))) / sum(from)
    w <- as.numeric(p_to / p_from)[as.integer(bin[from])]
    w[!is.finite(w)] <- 0
    w
}
n_eff <- function(w) sum(w)^2 / sum(w^2)

dm <- map(1:3, \(k) {
    c050 <- counts$cnt[[which(counts$tag == "mp050" & counts$k == k)]]
    c100 <- counts$cnt[[which(counts$tag == "mp100" & counts$k == k)]]
    si <- snp_info(c050)
    s100 <- snp_info(c100)
    dep_bin <- cut(si$rel_dep, unique(quantile(si$rel_dep, 0:20 / 20)),
                   include.lowest = TRUE)
    maf_bin <- cut(si$maf, unique(quantile(si$maf, 0:5 / 5)),
                   include.lowest = TRUE)
    dm_bin <- interaction(dep_bin, maf_bin, drop = FALSE)
    cmp <- si$complete
    inc <- !cmp
    all_ <- rep(TRUE, length(cmp))
    #' mp100 SNPs binned on the mp050 breaks
    dep100 <- cut(s100$rel_dep, unique(quantile(si$rel_dep, 0:20 / 20)),
                  include.lowest = TRUE)
    maf100 <- cut(s100$maf, unique(quantile(si$maf, 0:5 / 5)),
                  include.lowest = TRUE)
    dm100 <- interaction(dep100, maf100, drop = FALSE)
    w100 <- function(b050, b100) {
        p_to <- table(b050) / length(b050)
        p_from <- table(b100) / length(b100)
        w <- as.numeric(p_to / p_from)[as.integer(b100)]
        w[!is.finite(w)] <- 0
        w
    }
    rc <- c050$ref[cmp, ]
    ac <- c050$alt[cmp, ]
    w_d <- match_weights(dep_bin, cmp, inc)
    w_dm <- match_weights(dm_bin, cmp, inc)
    w100_d <- w100(dep_bin, dep100)
    w100_dm <- w100(dm_bin, dm100)
    list(
        n_eff = c(compl_depth = n_eff(w_d), compl_depth_maf = n_eff(w_dm),
                  mp100_depth = n_eff(w100_d), mp100_depth_maf = n_eff(w100_dm)),
        delta = list(
            incomplete = mom[[k]]$delta$mp050_incomplete,
            complete = mom[[k]]$delta$mp050_complete,
            complete_w_depth = pair_delta(rc, ac, c050$hap, w_d),
            complete_w_depth_maf = pair_delta(rc, ac, c050$hap, w_dm),
            mp050_all = mom[[k]]$delta$mp050_all,
            mp100 = mom[[k]]$delta$mp100,
            mp100_w_depth = pair_delta(c100$ref, c100$alt, c100$hap, w100_d),
            mp100_w_depth_maf = pair_delta(c100$ref, c100$alt, c100$hap,
                                           w100_dm)))
})

cat("\n== 5. Depth-matched comparison\n")
cat("Effective SNPs after weighting (sub1-3):\n")
print(sapply(dm, \(x) round(x$n_eff)))

dm_pairs <- tribble(
    ~a,           ~b,
    "incomplete", "complete",
    "incomplete", "complete_w_depth",
    "incomplete", "complete_w_depth_maf",
    "mp050_all",  "mp100",
    "mp050_all",  "mp100_w_depth",
    "mp050_all",  "mp100_w_depth_maf")
cat("\nD between sets, matched subsamples:\n")
dm_pairs |>
    mutate(d = map2(a, b, \(i, j) {
        map_dbl(1:3, \(k) dist_off(dm[[k]]$delta[[i]], dm[[k]]$delta[[j]]))
    }),
    d_min = map_dbl(d, min), d_max = map_dbl(d, max)) |>
    select(-d) |>
    print(width = Inf)
cat("\nSubsample noise (D between subsamples of the same set):\n")
tibble(set = names(dm[[1]]$delta)) |>
    mutate(d = map(set, \(s) {
        map_dbl(list(c(1, 2), c(1, 3), c(2, 3)), \(x) {
            dist_off(dm[[x[1]]]$delta[[s]], dm[[x[2]]]$delta[[s]])
        })
    }),
    d_min = map_dbl(d, min), d_max = map_dbl(d, max)) |>
    select(-d) |>
    print(width = Inf)

#' Sensitivity to the number of depth bins (depth only).
cat("\nD after depth matching, by number of depth bins (min-max over subsamples):\n")
crossing(n_bins = c(10, 20, 40), k = 1:3) |>
    mutate(d = map2(n_bins, k, \(nb, k) {
        c050 <- counts$cnt[[which(counts$tag == "mp050" & counts$k == k)]]
        c100 <- counts$cnt[[which(counts$tag == "mp100" & counts$k == k)]]
        si <- snp_info(c050)
        s100 <- snp_info(c100)
        brk <- unique(quantile(si$rel_dep, 0:nb / nb))
        b050 <- cut(si$rel_dep, brk, include.lowest = TRUE)
        b100 <- cut(s100$rel_dep, brk, include.lowest = TRUE)
        w_c <- match_weights(b050, si$complete, !si$complete)
        p_to <- table(b050) / length(b050)
        p_from <- table(b100) / length(b100)
        w_1 <- as.numeric(p_to / p_from)[as.integer(b100)]
        w_1[!is.finite(w_1)] <- 0
        cmp <- si$complete
        c(inc_vs_compl = dist_off(
              dm[[k]]$delta$incomplete,
              pair_delta(c050$ref[cmp, ], c050$alt[cmp, ], c050$hap, w_c)),
          all_vs_mp100 = dist_off(
              dm[[k]]$delta$mp050_all,
              pair_delta(c100$ref, c100$alt, c100$hap, w_1)))
    })) |>
    unnest_wider(d) |>
    group_by(n_bins) |>
    summarize(across(c(inc_vs_compl, all_vs_mp100),
                     list(min = min, max = max))) |>
    print(width = Inf)

sink()
#' Printed tibbles pad lines with spaces; strip them.
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)

#'
#' Results (2026-11-05; `_data/baypass-omega-diffs.txt`):
#'
#' 1. Not scale: mean log eigenvalue -0.03, mean-diagonal ratio 0.97-0.98,
#'    and FMD barely drops after rescaling (0.57-0.59 vs. 0.59-0.61).
#' 2. The largest and smallest eigenvalues (0.67, 1.4) carry most of FMD,
#'    and SN-1-95_S59 has 41-46% of the squared loading on both. It has the
#'    lowest depth (median 18) and the most missing sites (14%), and is the
#'    missing pool in 30% of SNPs with one missing pool.
#' 3. Pairwise divergence from reads matches BayPass's (r = 0.996).
#'    Imposing the missing pattern on complete SNPs doesn't change it
#'    (D 0.005, within subsample noise 0.007-0.011), and complete mp050
#'    SNPs match mp100 (D 0.005-0.006). Incomplete SNPs differ (D 0.08).
#'    So the difference is which SNPs are kept, not missing data itself
#'    (for this estimator; BayPass's own handling isn't tested here).
#' 4. Incomplete SNPs are at higher depth in the pools that have data
#'    (median relative depth 1.08-1.18 vs. 1.02). Divergence depends on
#'    depth in complete SNPs too: lower at relative depth > 1.2 (ratio
#'    0.79), higher at < 0.8 (1.26). Mid-depth incomplete SNPs still differ
#'    from complete ones (D 0.048 vs. 0.016 for mid-depth complete SNPs),
#'    so depth explains part of the difference, not all.
#' 5. Reweighting to match depth: all mp050 SNPs vs. mp100 drops from
#'    D 0.029-0.032 to 0.009-0.011, the same as subsample noise
#'    (0.008-0.010); incomplete vs. complete drops from 0.073-0.078 to
#'    0.027-0.030, still above noise (0.009-0.012). Matching MAF as well
#'    changes nothing; 10, 20, or 40 depth bins give the same result. So
#'    the mp050 vs. mp100 difference in divergence is explained by the
#'    depth of the SNPs each keeps (shown for read-based divergence, not
#'    BayPass's Omega itself).
#'
