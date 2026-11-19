#'
#' Why do pools differ so much from their neighbours in time? (Owner,
#' 2026-11-19, after the ABC prior-predictive check failed:
#' `_data/abc-results.txt`.) Differentiation between pools 1-2 generations
#' apart is far above any simulation and doesn't build up with time, i.e.
#' each pool has its own deviation. This measures that deviation per pool
#' and checks what it tracks.
#'
#' Per-pool excess e_j: from the pairwise divergence between pools (moment
#' estimate from reads with sampling variance for the pool's haploid size
#' and depth removed, as in `_R/baypass-pod-diagnose.R`; mp050 BayPass
#' inputs of the first rerun, all 3 subsamples), fit
#'   Δ_jk = e_j + e_k + b |t_j - t_k|
#' by non-negative least squares, with t in generations. e_j is the part of
#' pool j's divergence from every other pool that doesn't depend on time:
#' extra variance beyond the binomial sampling of 2 x adults genomes. If it
#' comes from a pool carrying fewer independent genomes than assumed, it
#' corresponds to an effective number of genomes n_eff with
#' 1 / n_eff = 1 / (2 x adults) + e_j (roughly). Also BayPass's Ω
#' diagonal (mp050, subsample 1), which absorbs the same excess.
#'
#' Then Spearman correlations of e_j with the pool's QC metrics
#' (`_data/pool-qc.csv`: missing data, duplication rates, GC, mapping
#' rates, SNPs and private SNPs), haploid pool size, median depth, and
#' sample year (older samples were stored longer), and trap abundance.
#' 17 pools (SN-1-95
#' excluded), so only strong relations can show.
#'
#' Writes `_data/pool-excess.txt` and `_plots/pool-excess.pdf`.
#'

library(tidyverse)

ARCHIVE_18POOLS <- TRUE
source("_R/data-dir.R")

out_file <- "_data/pool-excess.txt"
sink(out_file, split = TRUE)

pools18 <- read_lines("_bash/samples/time-pools-18.txt")
keep <- pools18 != "SN-1-95_S59"
pools <- pools18[keep]
J <- length(pools)



# =============================================================================*
# ---- data and pairwise divergence ----
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
ref <- g[, seq(1, ncol(g), 2)][, keep]
dep <- (g[, seq(1, ncol(g), 2)] + g[, seq(2, ncol(g), 2)])[, keep]
hap <- hap[keep]
rm(g)
untar(tany_dir("baypass-init", "baypass-init-mp050-sub1.tar.gz"),
      files = "baypass-init-mp050-sub1/tany-init-mp050-sub1_mat_omega.out",
      exdir = tmp)
om <- as.matrix(read.table(file.path(tmp, "baypass-init-mp050-sub1",
                                     "tany-init-mp050-sub1_mat_omega.out")))[keep, keep]
unlink(tmp, recursive = TRUE)

p_all <- ifelse(dep > 0, ref / pmax(dep, 1), NA_real_)
e_all <- rowMeans(p_all, na.rm = TRUE)
v_all <- e_all * (1 - e_all)
cc <- function(dp, h) 1 / pmax(dp, 1) + 1 / h - 1 / (pmax(dp, 1) * h)
delta <- matrix(NA_real_, J, J)
for (j in 1:(J - 1)) {
    for (k in (j + 1):J) {
        b <- dep[, j] > 0 & dep[, k] > 0
        delta[j, k] <- delta[k, j] <-
            sum((p_all[b, j] - p_all[b, k])^2 -
                    v_all[b] * (cc(dep[b, j], hap[j]) + cc(dep[b, k], hap[k]))) /
            sum(v_all[b])
    }
}

info <- tibble(pool = pools) |>
    mutate(yy = as.integer(str_match(pool, "^SN-[12]-(\\d{2})")[, 2]),
           year = ifelse(yy < 70, 2000L, 1900L) + yy,
           t = 2L * (year - 1977L) +
               as.integer(str_match(pool, "^SN-([12])-")[, 2]) - 1L,
           haploid = hap,
           median_depth = apply(ifelse(dep > 0, dep, NA), 2, median,
                                na.rm = TRUE))



# =============================================================================*
# ---- fit ----
# =============================================================================*


pr <- t(combn(J, 2))
X <- matrix(0, nrow(pr), J + 1)
X[cbind(seq_len(nrow(pr)), pr[, 1])] <- 1
X[cbind(seq_len(nrow(pr)), pr[, 2])] <- 1
X[, J + 1] <- abs(info$t[pr[, 1]] - info$t[pr[, 2]])
y <- delta[pr]
fit <- nnls::nnls(X, y)
e <- fit$x[1:J]
b_time <- fit$x[J + 1]
r2 <- 1 - sum((y - X %*% fit$x)^2) / sum((y - mean(y))^2)
#' Without per-pool terms, for comparison:
r2_time <- summary(lm(y ~ X[, J + 1]))$r.squared

info <- info |>
    mutate(excess = e,
           omega_diag = diag(om),
           n_eff = 1 / (1 / haploid + excess))
cat("== Pairwise divergence = e_j + e_k + b x generations apart\n")
cat("b (per generation):", signif(b_time, 3), "; R^2", signif(r2, 4),
    "(time alone:", signif(r2_time, 3), "); residual SD",
    signif(sd(y - X %*% fit$x), 3), "vs. SD of the divergences",
    signif(sd(y), 3), "\n")
cat("Implied effective size from b alone: Ne ~ 1 / (2 b) =",
    signif(1 / (2 * b_time), 3), "\n\n")
info |>
    arrange(desc(excess)) |>
    mutate(across(c(excess, omega_diag), \(x) signif(x, 3)),
           n_eff = round(n_eff, 1)) |>
    select(pool, year, haploid, median_depth, excess, n_eff, omega_diag) |>
    as.data.frame() |>
    print()
cat("\nSpearman correlation of excess with the BayPass Omega diagonal:",
    round(cor(info$excess, info$omega_diag, method = "spearman"), 3), "\n")



# =============================================================================*
# ---- what does it track? ----
# =============================================================================*


qc <- read_csv("_data/pool-qc.csv", show_col_types = FALSE) |>
    filter(set == "time") |>
    select(pool, prop_missing, picard_dup, fastp_dup, fastp_gc, target_rate,
           mitolike_rate, map_rate, n_snps, n_private)
#' Trap abundance: log(1 + N) in the pool's generation and hm_logN_6 (the
#' BayPass covariate; `_data/archive-18pools/baypass-covariates.csv`). If
#' the excess comes from sampling related individuals, it should be larger
#' when there are fewer adults to sample from.
abund <- read_csv("_data/archive-18pools/baypass-covariates.csv",
                  show_col_types = FALSE) |>
    select(pool, logN, hm_logN_6)
dd <- info |> left_join(qc, by = "pool") |> left_join(abund, by = "pool")
stopifnot(nrow(dd) == J, !anyNA(dd$picard_dup))
vars <- c("prop_missing", "picard_dup", "fastp_dup", "fastp_gc", "target_rate",
          "mitolike_rate", "map_rate", "n_snps", "n_private", "haploid",
          "median_depth", "year", "logN", "hm_logN_6")
cors <- map(vars, \(v) {
    ct <- suppressWarnings(cor.test(dd$excess, dd[[v]], method = "spearman"))
    tibble(variable = v, rho = ct$estimate, p = ct$p.value)
}) |>
    list_rbind()
cat("\n== Spearman correlation of the excess with pool variables (17 pools)\n")
print(as.data.frame(cors |> mutate(rho = round(rho, 3), p = signif(p, 2)) |>
                        arrange(p)))

p1 <- dd |>
    select(pool, excess, all_of(vars)) |>
    pivot_longer(all_of(vars)) |>
    ggplot(aes(value, excess)) +
    geom_point() +
    facet_wrap(~ name, scales = "free_x") +
    labs(x = NULL, y = "Per-pool excess divergence") +
    theme_bw()
p2 <- tibble(gap = X[, J + 1], delta = y,
             pair_excess = X[, 1:J] %*% e) |>
    ggplot(aes(gap, delta)) +
    geom_point(aes(color = pair_excess)) +
    geom_abline(intercept = 0, slope = b_time, linetype = 2) +
    scale_color_viridis_c() +
    labs(x = "Generations apart", y = "Pairwise divergence (reads)",
         color = "e_j + e_k") +
    theme_classic()
pdf("_plots/pool-excess.pdf", width = 9, height = 7)
print(p2)
print(p1)
invisible(dev.off())

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
