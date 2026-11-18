#'
#' Power test of candidate statistics on PODs with spiked-in associations
#' (owner, 2026-11-18), before any real Bayes factor is read. Compares
#'   bf      median BF (dB) over 10 runs (the rule fixed 2026-11-07)
#'   z       median |M_Beta| / SD_Beta over 10 runs
#'   bf_z    median BF, ties broken by median z
#' by how many spiked SNPs each calls, and how many null SNPs, at Storey
#' q <= 0.05 and q <= 0.10, with empirical p-values against the existing
#' null POD of the same setting and covariate
#' (`tany_dir("baypass-pod", "pod-null-stats.rds")`, from
#' `_R/baypass-pod-stat-compare.R`).
#'
#' Spiked POD: setting mp075 (subsample-1 Ω, π, depths, pool sizes, and
#' filters, as in `_R/baypass-pod-inputs.R`), covariate `hm_logN_6`.
#' 100,000 SNPs after filters: 99,000 null and 1,000 associated. As in
#' `simulate.baypass()` (`_R/baypass-utils.R`), an associated SNP's
#' population frequencies get beta x c_j added before truncation, with c
#' the covariate scaled to mean 0 and SD 1 (BayPass scales covariates), and
#' |beta| log-uniform on 0.01-0.2 (frequency change per SD of the
#' covariate; weak to strong), random sign. Run through the covariate model
#' as the real runs (`03-select-time/06-baypass-pod-spike`, 10 runs with
#' the real runs' seeds).
#'
#' Part 1 (inputs): writes `tany_dir("baypass-pod",
#' "baypass-pod-spike-inputs.tar.gz")` (upload to
#' `/staging/l/lnell/dna/baypass/`). Part 2 (results; when the 10 outputs
#' are in `tany_dir("baypass-pod")`). Both report to
#' `_data/baypass-pod-spike.txt`.
#'

library(tidyverse)

source("_R/data-dir.R")

local({
    for (x in parse("_R/baypass-pod-inputs.R")) {
        if (is.call(x) && identical(x[[1]], as.name("<-")) &&
            is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function"))) {
            eval(x, envir = globalenv())
        }
    }
})

out_file <- "_data/baypass-pod-spike.txt"
sink(out_file, split = TRUE)

tag <- "mp075"
cov_name <- "hm_logN_6"
n_null <- 99000L
n_spike <- 1000L
beta_range <- c(0.01, 0.2)
n_runs <- 10
spike_tar <- tany_dir("baypass-pod", "baypass-pod-spike-inputs.tar.gz")

set.seed(1880354617)



# =============================================================================*
# ---- part 1: inputs ----
# =============================================================================*


if (!file.exists(spike_tar)) {
    src <- get_pod_src(tag)
    J <- length(src$hap)
    cv <- scan(file.path("_bash/03-select-time/covariates",
                         paste0(cov_name, ".txt")), quiet = TRUE)
    stopifnot(length(cv) == J)
    cs <- (cv - mean(cv)) / sd(cv)

    #' As `simulate_snps()`, plus beta x c for associated SNPs.
    sim_assoc <- function(n, beta) {
        i_pi <- sample.int(length(src$pi_hat), n, replace = TRUE)
        i_dep <- sample.int(nrow(src$dep), n, replace = TRUE)
        pi_ <- src$pi_hat[i_pi]
        dep <- src$dep[i_dep, ]
        z <- matrix(rnorm(n * J), J, n)
        alpha <- t(t(chol(src$omega)) %*% z) * sqrt(pi_ * (1 - pi_)) + pi_ +
            outer(beta, cs)
        alpha <- pmin(pmax(alpha, 0), 1)
        hap <- matrix(src$hap, n, J, byrow = TRUE)
        y <- matrix(rbinom(n * J, hap, alpha), n, J)
        ref <- matrix(rbinom(n * J, dep, y / hap), n, J)
        list(ref = ref, dep = dep, beta = beta)
    }
    draw_until <- function(n, beta_fun) {
        ref <- dep <- NULL
        beta <- numeric(0)
        while (length(beta) < n) {
            b <- beta_fun(2 * n)
            s <- sim_assoc(2 * n, b)
            k <- passes_filters(s$ref, s$dep)
            ref <- rbind(ref, s$ref[k, ])
            dep <- rbind(dep, s$dep[k, ])
            beta <- c(beta, b[k])
        }
        list(ref = ref[1:n, ], dep = dep[1:n, ], beta = beta[1:n])
    }
    nul <- draw_until(n_null, \(n) rep(0, n))
    spk <- draw_until(n_spike, \(n) {
        exp(runif(n, log(beta_range[1]), log(beta_range[2]))) *
            sample(c(-1, 1), n, replace = TRUE)
    })
    #' Shuffle so associated SNPs aren't together:
    o <- sample.int(n_null + n_spike)
    ref <- rbind(nul$ref, spk$ref)[o, ]
    dep <- rbind(nul$dep, spk$dep)[o, ]
    truth <- tibble(snp = seq_along(o),
                    beta = c(nul$beta, spk$beta)[o])

    out_dir <- file.path(tempdir(), "baypass-pod-spike-inputs")
    d <- file.path(out_dir, "baypass-podspike")
    dir.create(d, recursive = TRUE)
    g <- matrix(0L, nrow(ref), 2 * J)
    g[, seq(1, ncol(g), 2)] <- ref
    g[, seq(2, ncol(g), 2)] <- dep - ref
    data.table::fwrite(as.data.frame(g), file.path(d, "tany-pod.genobaypass"),
                       sep = " ", col.names = FALSE)
    write_lines(paste(src$hap, collapse = " "), file.path(d, "tany-pod.poolsize"))
    write_csv(truth, file.path(d, "truth.csv"))
    write_rds(truth, tany_dir("baypass-pod", "pod-spike-truth.rds"),
              compress = "gz")
    old <- setwd(out_dir)
    tar(spike_tar, files = "baypass-podspike", compression = "gzip",
        tar = "internal")
    setwd(old)
    cat("== Part 1: spiked POD (", tag, ", ", cov_name, "): ", n_null,
        " null + ", n_spike, " associated SNPs, |beta| log-uniform ",
        beta_range[1], "-", beta_range[2], "\n", sep = "")
    cat("Wrote", spike_tar, "\n")
}



# =============================================================================*
# ---- part 2: results ----
# =============================================================================*


res_name <- function(r) sprintf("baypass-podspike-run%i", r)
have <- file.exists(tany_dir("baypass-pod", paste0(res_name(seq_len(n_runs)),
                                                   ".tar")))
if (all(have)) {
    truth <- read_rds(tany_dir("baypass-pod", "pod-spike-truth.rds"))
    runs <- map(seq_len(n_runs), \(r) {
        nm <- res_name(r)
        member <- file.path(nm, sprintf("tany-podspike-run%i_summary_betai.out.gz", r))
        tmp <- tempfile()
        untar(tany_dir("baypass-pod", paste0(nm, ".tar")), files = member,
              exdir = tmp)
        x <- data.table::fread(cmd = paste0("gzip -dc '", file.path(tmp, member),
                                            "'"),
                               select = c("M_Beta", "SD_Beta", "BF(dB)"))
        cbind(bf = x[["BF(dB)"]],
              z = ifelse(x$SD_Beta > 0, abs(x$M_Beta) / x$SD_Beta, 0))
    })
    obs <- tibble(snp = truth$snp, beta = truth$beta,
                  bf = apply(sapply(runs, \(x) x[, "bf"]), 1, median),
                  z = apply(sapply(runs, \(x) x[, "z"]), 1, median))
    nul <- read_rds(tany_dir("baypass-pod", "pod-null-stats.rds")) |>
        filter(tag == !!tag, covariate == cov_name)
    n0 <- nrow(nul)

    #' Empirical p: (1 + null SNPs at least as extreme) / (n0 + 1).
    p_one <- function(x, null) {
        (1 + n0 - findInterval(x, sort(null), left.open = TRUE)) / (n0 + 1)
    }
    #' BF with ties broken by z, for SNPs at or above the null's largest
    #' BF (the ceiling, where the ties are): a null SNP is at least as
    #' extreme if its BF is larger, or equal and its z >= the SNP's.
    p_tie <- function(bf, z) {
        p <- p_one(bf, nul$bf)
        top <- which(bf >= max(nul$bf))
        for (i in top) {
            ge <- sum(nul$bf > bf[i] | (nul$bf == bf[i] & nul$z >= z[i]))
            p[i] <- (1 + ge) / (n0 + 1)
        }
        p
    }
    obs <- obs |>
        mutate(p_bf = p_one(bf, nul$bf),
               p_z = p_one(z, nul$z),
               p_bf_z = p_tie(bf, z))
    obs <- obs |>
        mutate(q_bf = qvalue::qvalue(p_bf)$qvalues,
               q_z = qvalue::qvalue(p_z)$qvalues,
               q_bf_z = qvalue::qvalue(p_bf_z)$qvalues)

    cat("\n== Part 2: calls at q <= 0.05 and q <= 0.10 (", n_spike,
        " associated, ", n_null, " null SNPs)\n", sep = "")
    calls <- crossing(stat = c("bf", "z", "bf_z"), q_max = c(0.05, 0.10)) |>
        mutate(x = map2(stat, q_max, \(s, qm) {
            called <- obs[[paste0("q_", s)]] <= qm
            tibble(true_pos = sum(called & obs$beta != 0),
                   false_pos = sum(called & obs$beta == 0),
                   fdr = sum(called & obs$beta == 0) / max(sum(called), 1),
                   power = sum(called & obs$beta != 0) / n_spike)
        })) |>
        unnest(x)
    calls |>
        mutate(across(c(fdr, power), \(x) round(x, 3))) |>
        as.data.frame() |>
        print()

    cat("\nPower at q <= 0.10 by |beta| (quintile bins of the spiked SNPs):\n")
    obs |>
        filter(beta != 0) |>
        mutate(bin = cut(abs(beta), quantile(abs(beta), 0:5 / 5),
                         include.lowest = TRUE)) |>
        group_by(bin) |>
        summarize(n = n(), bf = mean(q_bf <= 0.1), z = mean(q_z <= 0.1),
                  bf_z = mean(q_bf_z <= 0.1)) |>
        mutate(across(where(is.double), \(x) round(x, 3))) |>
        as.data.frame() |>
        print()
    cat("\nSpiked SNPs at the BF ceiling (median):",
        sum(obs$beta != 0 & obs$bf >= 59.95), "; null SNPs:",
        sum(obs$beta == 0 & obs$bf >= 59.95), "\n")
    write_rds(obs, tany_dir("baypass-pod", "pod-spike-results.rds"),
              compress = "gz")
} else {
    cat("\nPart 2: waiting for", sum(!have), "of", length(have), "outputs\n")
}

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
