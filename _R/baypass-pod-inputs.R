#'
#' Pseudo-observed data (PODs) to calibrate the BayPass Bayes factor (BF)
#' cutoff (BayPass manual v3.1, section 4.1; Gautier 2015). For each
#' missing-data setting, SNPs are simulated under the core model, with no
#' association with any covariate, from that setting's estimates in
#' `01-baypass-init` (subsample 1):
#'   - Omega;
#'   - ancestral allele frequencies pi: each simulated SNP takes the
#'     posterior mean of pi (`M_P`) of a random real SNP;
#'   - haploid pool sizes and read depths, each simulated SNP taking the
#'     depths (including 0 = masked pool) of a random real SNP.
#' The model is that of `simulate.baypass()` (`_R/baypass-utils.R`), which
#' draws pi from the fitted Beta(a_pi, b_pi) instead. Resampling the real
#' estimates keeps the frequency spectrum of the filtered SNPs, which a
#' two-parameter Beta need not match; the frequency checks below compare
#' PODs with the real SNPs.
#' Simulated SNPs then pass the same filters as the real ones in
#' `00-baypass-inputs.sh`: pooled minor-allele read frequency >= 0.01
#' (poolfstat default), and mean reference read proportion over pools with
#' data in (0.1, 0.9). (The real SNPs were also called polymorphic by SNAPE
#' in at least one pool; that step isn't mimicked.)
#'
#' `n_pod` = 200,000 SNPs per setting: the real settings test 0.53-1.02
#' million SNPs, so a null tail probability of 1e-4 (~100 false positives
#' among 1 million SNPs) is estimated from ~20 POD SNPs. At the measured
#' ~24 h per 341,286 SNPs for the covariate model, a POD job takes ~14 h.
#'
#' Reads `tany_dir("baypass-init")` and `tany_dir("baypass-inputs")`.
#' Writes `tany_dir("baypass-pod", "baypass-pod-inputs.tar.gz")` (upload to
#' `/staging/l/lnell/dna/baypass/` for `03-select-time/04-baypass-pod`) and
#' `_data/baypass-pod-inputs.txt` (POD vs. real data checks).
#'

library(tidyverse)

source("_R/data-dir.R")

out_file <- "_data/baypass-pod-inputs.txt"
sink(out_file, split = TRUE)

mp_tags <- c("mp050", "mp075", "mp100")
n_pod <- 200000L
#' Simulated per batch, before filtering (most pass):
n_batch <- 100000L

set.seed(1771348075)



# =============================================================================*
# ---- inputs ----
# =============================================================================*


#' Omega, posterior means of pi, Beta(a_pi, b_pi) parameters (reported
#' only), read depths (SNPs x pools), and haploid pool sizes for one
#' setting (subsample 1).
get_pod_src <- function(tag) {
    tmp <- tempfile()
    nm <- sprintf("baypass-init-%s-sub1", tag)
    pre <- sprintf("%s/tany-init-%s-sub1", nm, tag)
    untar(tany_dir("baypass-init", paste0(nm, ".tar.gz")),
          files = paste0(pre, c("_mat_omega.out", "_summary_beta_params.out",
                                "_summary_pi_xtx.out")),
          exdir = tmp)
    omega <- file.path(tmp, paste0(pre, "_mat_omega.out")) |>
        read_table(col_names = FALSE, col_types = cols()) |>
        as.matrix() |>
        unname()
    bp <- file.path(tmp, paste0(pre, "_summary_beta_params.out")) |>
        read_table(col_types = cols())
    pi_hat <- file.path(tmp, paste0(pre, "_summary_pi_xtx.out")) |>
        data.table::fread(data.table = FALSE) |>
        getElement("M_P")
    d <- sprintf("baypass-inputs_%s", tag)
    fns <- file.path(d, c("tany.genobaypass.sub1", "tany.poolsize"))
    untar(tany_dir("baypass-inputs", sprintf("baypass-inputs_%s_3.tar.gz", tag)),
          files = fns, exdir = tmp)
    g <- data.table::fread(file.path(tmp, fns[1]), header = FALSE) |>
        as.matrix() |>
        unname()
    hap <- read_lines(file.path(tmp, fns[2])) |>
        str_split_1(" ") |>
        as.integer()
    stopifnot(length(pi_hat) == nrow(g))
    list(omega = omega,
         pi_hat = pi_hat,
         beta_pi = bp$Mean[match(c("a_beta_pi", "b_beta_pi"), bp$PARAM)],
         ref = g[, seq(1, ncol(g), 2)],
         dep = g[, seq(1, ncol(g), 2)] + g[, seq(2, ncol(g), 2)],
         hap = hap)
}

#' Same filters as `00-baypass-inputs.sh` (see header). `ref`, `dep`:
#' SNPs x pools read counts.
passes_filters <- function(ref, dep) {
    has <- dep > 0
    pooled <- rowSums(ref) / rowSums(dep)
    mean_ref <- rowSums(ifelse(has, ref / pmax(dep, 1), 0)) / rowSums(has)
    pmin(pooled, 1 - pooled) >= 0.01 & mean_ref > 0.1 & mean_ref < 0.9
}



# =============================================================================*
# ---- simulate ----
# =============================================================================*


#' Core model, as in `simulate.baypass()`: population frequencies
#' alpha_j = pi + sqrt(pi (1 - pi)) (L z)_j, with L L' = Omega (Cholesky)
#' and z standard normal, truncated to [0, 1]; allele counts in the pool
#' ~ Binomial(haploid size, alpha_j); reference reads ~ Binomial(depth,
#' count / haploid size).
simulate_snps <- function(src, n) {
    J <- length(src$hap)
    i_pi <- sample.int(length(src$pi_hat), n, replace = TRUE)
    i_dep <- sample.int(nrow(src$dep), n, replace = TRUE)
    pi_ <- src$pi_hat[i_pi]
    dep <- src$dep[i_dep, ]
    z <- matrix(rnorm(n * J), J, n)
    alpha <- t(t(chol(src$omega)) %*% z) * sqrt(pi_ * (1 - pi_)) + pi_
    alpha <- pmin(pmax(alpha, 0), 1)
    hap <- matrix(src$hap, n, J, byrow = TRUE)
    y <- matrix(rbinom(n * J, hap, alpha), n, J)
    ref <- matrix(rbinom(n * J, dep, y / hap), n, J)
    list(ref = ref, dep = dep)
}

simulate_pod <- function(src) {
    ref <- NULL
    dep <- NULL
    n_sim <- 0L
    while (is.null(ref) || nrow(ref) < n_pod) {
        sim <- simulate_snps(src, n_batch)
        n_sim <- n_sim + n_batch
        keep <- passes_filters(sim$ref, sim$dep)
        ref <- rbind(ref, sim$ref[keep, ])
        dep <- rbind(dep, sim$dep[keep, ])
    }
    list(ref = ref[1:n_pod, ], dep = dep[1:n_pod, ], n_sim = n_sim,
         n_pass = nrow(ref))
}

pods <- mp_tags |>
    set_names() |>
    map(\(tag) {
        src <- get_pod_src(tag)
        pod <- simulate_pod(src)
        list(src = src, pod = pod)
    })



# =============================================================================*
# ---- checks: PODs vs. real data ----
# =============================================================================*


#' Summary of a set of SNPs (ref, dep: SNPs x pools).
snp_summ <- function(ref, dep) {
    has <- dep > 0
    mean_ref <- rowSums(ifelse(has, ref / pmax(dep, 1), 0)) / rowSums(has)
    tibble(snps = nrow(ref),
           prop_missing = mean(!has),
           median_depth = median(dep[has]),
           mean_ref_q10 = quantile(mean_ref, 0.1),
           mean_ref_q50 = quantile(mean_ref, 0.5),
           mean_ref_q90 = quantile(mean_ref, 0.9),
           maf_q50 = median(pmin(mean_ref, 1 - mean_ref)))
}

cat("== Simulation parameters (from 01-baypass-init, subsample 1)\n")
imap(pods, \(x, tag) {
    tibble(tag = tag, a_beta_pi = x$src$beta_pi[1], b_beta_pi = x$src$beta_pi[2],
           simulated = x$pod$n_sim, passed_filters = x$pod$n_pass,
           kept = nrow(x$pod$ref))
}) |>
    list_rbind() |>
    as.data.frame() |>
    print()

#' The real SNPs passed these filters already; this checks the function.
cat("\nReal SNPs (subsample 1) passing the filters (should be 1):\n")
print(map_dbl(pods, \(x) mean(passes_filters(x$src$ref, x$src$dep))))

cat("\n== Real (subsample 1) vs. POD SNPs\n")
imap(pods, \(x, tag) {
    bind_rows(snp_summ(x$src$ref, x$src$dep) |> mutate(set = "real"),
              snp_summ(x$pod$ref, x$pod$dep) |> mutate(set = "pod")) |>
        mutate(tag = tag)
}) |>
    list_rbind() |>
    select(tag, set, everything()) |>
    mutate(across(where(is.double), \(x) signif(x, 3))) |>
    as.data.frame() |>
    print()



# =============================================================================*
# ---- write ----
# =============================================================================*


#' One folder per setting, as `00-baypass-inputs` writes them: genotype
#' file (ref and alt read counts per pool) and pool sizes.
out_dir <- file.path(tempdir(), "baypass-pod-inputs")
dir.create(out_dir)
iwalk(pods, \(x, tag) {
    d <- file.path(out_dir, paste0("baypass-pod_", tag))
    dir.create(d)
    g <- matrix(0L, n_pod, 2 * ncol(x$pod$ref))
    g[, seq(1, ncol(g), 2)] <- x$pod$ref
    g[, seq(2, ncol(g), 2)] <- x$pod$dep - x$pod$ref
    data.table::fwrite(as.data.frame(g), file.path(d, "tany-pod.genobaypass"),
                       sep = " ", col.names = FALSE)
    write_lines(paste(x$src$hap, collapse = " "),
                file.path(d, "tany-pod.poolsize"))
})

dir.create(tany_dir("baypass-pod"), showWarnings = FALSE)
tar_file <- tany_dir("baypass-pod", "baypass-pod-inputs.tar.gz")
old <- setwd(out_dir)
tar(tar_file, files = paste0("baypass-pod_", mp_tags), compression = "gzip",
    tar = "internal")
setwd(old)
cat("\nWrote", tar_file, "\n")

sink()
#' Printed tables pad lines with spaces; strip them.
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
