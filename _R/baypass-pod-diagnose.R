#'
#' Why the Ω re-estimated from the PODs differs from the Ω they were
#' simulated from (FMD 1.04-1.13, diagonals 0.83-1.05x;
#' `_data/baypass-pod-threshold.txt`, step 1). Simulation-only check, no
#' BayPass runs (owner, 2026-11-17): does the POD simulation itself keep
#' the drift variance that Ω implies, and if not, which step loses it?
#'
#' Setting mp075, subsample 1 (the source of the mp075 POD). Measure:
#' moment-based pairwise divergence Δ_jl = Ω_jj + Ω_ll - 2 Ω_jl, from read
#' counts with sampling variance removed (`pair_delta`, as in
#' `_R/baypass-omega-differences.R`), compared with Δ computed directly
#' from the source Ω. Reported per variant: mean over pool pairs of
#' Δ(variant) / Δ(Ω), its range over pools (each pool's mean over its
#' pairs), and the correlation over pairs. The real reads give the
#' reference: what the moment method returns for data whose BayPass Ω is
#' the source Ω.
#'
#' Variants (200,000 simulated SNPs each, before filters):
#'   pod        the POD method (`_R/baypass-pod-inputs.R`): π resampled
#'              from the posterior means (M_P), population frequencies
#'              truncated to [0, 1], real depths, then the filters
#'   no_filter  the same, before the filters
#'   beta_pi    π from the fitted Beta(a_pi, b_pi) (as `simulate.baypass`)
#'   alpha      population frequencies themselves (before pool and read
#'              sampling), with and without truncation: Δ from α needs
#'              no sampling correction
#'
#' Writes `_data/baypass-pod-diagnose.txt`.
#'

library(tidyverse)

source("_R/data-dir.R")

#' Functions from `_R/baypass-pod-inputs.R` (only its function
#' definitions are evaluated):
local({
    for (x in parse("_R/baypass-pod-inputs.R")) {
        if (is.call(x) && identical(x[[1]], as.name("<-")) &&
            is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function"))) {
            eval(x, envir = globalenv())
        }
    }
})

out_file <- "_data/baypass-pod-diagnose.txt"
sink(out_file, split = TRUE)

set.seed(1720664413)
tag <- "mp075"
n_sim <- 200000L

src <- get_pod_src(tag)
J <- length(src$hap)
cat("Setting", tag, "subsample 1:", nrow(src$ref), "real SNPs,", J,
    "pools\n")
cat("Beta(a_pi, b_pi):", paste(signif(src$beta_pi, 4), collapse = ", "), "\n")



# =============================================================================*
# ---- divergence ----
# =============================================================================*


#' As in `_R/baypass-omega-differences.R` (section 3).
pair_delta <- function(ref, alt, hap) {
    dep <- ref + alt
    has <- dep > 0
    p <- ifelse(has, ref / pmax(dep, 1), NA_real_)
    e <- rowMeans(p, na.rm = TRUE)
    v <- e * (1 - e)
    cc <- sweep(1 / pmax(dep, 1), 2, 1 / hap, "+") -
        sweep(1 / pmax(dep, 1), 2, hap, "/")
    out <- matrix(0, J, J)
    for (j in 1:(J - 1)) {
        for (l in (j + 1):J) {
            b <- has[, j] & has[, l]
            out[j, l] <- out[l, j] <-
                sum((p[b, j] - p[b, l])^2 - v[b] * (cc[b, j] + cc[b, l])) /
                sum(v[b])
        }
    }
    out
}
#' Δ from population frequencies (no sampling), centered on π:
#' sum_i (a_ij - a_il)^2 / sum_i π_i (1 - π_i).
alpha_delta <- function(alpha, pi_) {
    out <- matrix(0, J, J)
    v <- sum(pi_ * (1 - pi_))
    for (j in 1:(J - 1)) {
        for (l in (j + 1):J) {
            out[j, l] <- out[l, j] <- sum((alpha[, j] - alpha[, l])^2) / v
        }
    }
    out
}
delta_of_omega <- function(om) outer(diag(om), diag(om), "+") - 2 * om

d_om <- delta_of_omega(src$omega)
off <- upper.tri(d_om)
compare <- function(d, label) {
    r <- d / d_om
    diag(r) <- NA
    per_pool <- rowMeans(r, na.rm = TRUE)
    tibble(variant = label, mean_ratio = mean(r[off]),
           pool_min = min(per_pool), pool_max = max(per_pool),
           cor = cor(d[off], d_om[off]))
}



# =============================================================================*
# ---- variants ----
# =============================================================================*


#' Simulation as `simulate_snps()`, with the π source and truncation as
#' options, returning population frequencies too.
sim_variant <- function(n, pi_from = c("mp", "beta"), truncate = TRUE) {
    pi_from <- match.arg(pi_from)
    pi_ <- if (pi_from == "mp") {
        src$pi_hat[sample.int(length(src$pi_hat), n, replace = TRUE)]
    } else {
        rbeta(n, src$beta_pi[1], src$beta_pi[2])
    }
    dep <- src$dep[sample.int(nrow(src$dep), n, replace = TRUE), ]
    z <- matrix(rnorm(n * J), J, n)
    alpha_raw <- t(t(chol(src$omega)) %*% z) * sqrt(pi_ * (1 - pi_)) + pi_
    alpha <- pmin(pmax(alpha_raw, 0), 1)
    a_use <- if (truncate) alpha else alpha_raw
    hap <- matrix(src$hap, n, J, byrow = TRUE)
    y <- matrix(rbinom(n * J, hap, alpha), n, J)
    ref <- matrix(rbinom(n * J, dep, y / hap), n, J)
    list(ref = ref, dep = dep, pi = pi_, alpha = a_use,
         truncated = mean(alpha_raw < 0 | alpha_raw > 1))
}

res <- list()
res$real <- compare(pair_delta(src$ref, src$dep - src$ref, src$hap),
                    "real reads (sub1)")

s_mp <- sim_variant(n_sim, "mp")
keep <- passes_filters(s_mp$ref, s_mp$dep)
res$pod <- compare(pair_delta(s_mp$ref[keep, ], (s_mp$dep - s_mp$ref)[keep, ],
                              src$hap), "pod (M_P, truncated, filtered)")
res$no_filter <- compare(pair_delta(s_mp$ref, s_mp$dep - s_mp$ref, src$hap),
                         "pod before filters")
res$alpha_trunc <- compare(alpha_delta(s_mp$alpha, s_mp$pi),
                           "alpha, M_P, truncated")
s_raw <- sim_variant(n_sim, "mp", truncate = FALSE)
res$alpha_raw <- compare(alpha_delta(s_raw$alpha, s_raw$pi),
                         "alpha, M_P, not truncated")

s_b <- sim_variant(n_sim, "beta")
keep_b <- passes_filters(s_b$ref, s_b$dep)
res$beta <- compare(pair_delta(s_b$ref[keep_b, ], (s_b$dep - s_b$ref)[keep_b, ],
                               src$hap), "pod with Beta pi (filtered)")
res$alpha_beta <- compare(alpha_delta(s_b$alpha, s_b$pi),
                          "alpha, Beta pi, truncated")

cat("\n== Δ(variant) / Δ(source Ω), over pool pairs\n")
cat("(mean_ratio: mean over pairs; pool_min/max: range of each pool's mean",
    "over its pairs; cor: correlation over pairs)\n")
bind_rows(res) |>
    mutate(across(where(is.double), \(x) round(x, 3))) |>
    as.data.frame() |>
    print()

cat("\nShare of population frequencies truncated to [0, 1]: M_P pi ",
    signif(s_mp$truncated, 3), "; Beta pi ", signif(s_b$truncated, 3), "\n",
    sep = "")
cat("Share of simulated SNPs passing filters: M_P pi ", signif(mean(keep), 3),
    "; Beta pi ", signif(mean(keep_b), 3), "\n", sep = "")
cat("pi spread, sd: real M_P ", signif(sd(src$pi_hat), 3), "; Beta ",
    signif(sd(s_b$pi), 3), "\n", sep = "")

#' Per pool: each pool's mean Δ ratio over its pairs (real reads, POD
#' reads), next to the ratio of BayPass's Ω diagonal re-estimated from the
#' full mp075 POD to the source Ω (step 1 of
#' `_R/baypass-pod-threshold.R`). Pool order is that of
#' `_bash/samples/time-pools-18.txt` (as in the BayPass inputs).
pool_ratio <- function(d) {
    r <- d / d_om
    diag(r) <- NA
    rowMeans(r, na.rm = TRUE)
}
nm <- sprintf("baypass-pod-%s-core-run1", tag)
tmp <- tempfile()
untar(tany_dir("baypass-pod", paste0(nm, ".tar")),
      files = file.path(nm, sprintf("tany-pod-%s-core-run1_mat_omega.out", tag)),
      exdir = tmp)
om_re <- as.matrix(read.table(file.path(tmp, nm, sprintf(
    "tany-pod-%s-core-run1_mat_omega.out", tag))))
dep_real <- src$dep
dep_real[dep_real == 0] <- NA
cat("\n== Per pool\n")
pp <- tibble(pool = read_lines("_bash/samples/time-pools-18.txt"),
       median_depth = apply(dep_real, 2, median, na.rm = TRUE),
       haploid = src$hap,
       omega_diag = diag(src$omega),
       real_delta_ratio = pool_ratio(pair_delta(src$ref, src$dep - src$ref,
                                                src$hap)),
       pod_delta_ratio = pool_ratio(pair_delta(
           s_mp$ref[keep, ], (s_mp$dep - s_mp$ref)[keep, ], src$hap)),
       baypass_pod_diag_ratio = diag(om_re) / diag(src$omega))
pp |>
    mutate(across(where(is.double), \(x) signif(x, 3))) |>
    as.data.frame() |>
    print()
cat("\nSpearman correlation, source Omega diagonal vs. BayPass re-estimated",
    "/ source diagonal:",
    signif(cor(pp$omega_diag, pp$baypass_pod_diag_ratio, method = "spearman"),
           3), "\n")

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
