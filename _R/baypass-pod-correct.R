#'
#' Test of a bias correction for the POD Ω (owner, 2026-11-17), at 100,000
#' SNPs per setting before any full-size rerun.
#'
#' Problem (`_data/baypass-pod-threshold.txt` step 1;
#' `_data/baypass-pod-diagnose.txt`): BayPass's Ω re-estimated from the
#' PODs is 6-16% low on the diagonal for low-drift pools, although the
#' simulation keeps the divergence the source Ω implies. Correction: one
#' step, simulate from Ω* = D Ω D, with D diagonal, D_jj = sqrt(Ω_jj /
#' Ω_re,jj), Ω the subsample-1 Ω (source), and Ω_re BayPass's re-estimate
#' from the full POD (core model, `04-baypass-pod`). This rescales each
#' pool's drift variance by its shortfall and keeps the correlations among
#' pools (which BayPass recovered: off-diagonal correlation 0.997-0.998),
#' and Ω* stays positive definite. The additive step Ω + (Ω - Ω_re) was
#' tried first and isn't positive definite for mp050 (smallest eigenvalue
#' -0.00019). If BayPass's shrinkage is about the same near Ω*,
#' re-estimating from a POD simulated from Ω* returns ~Ω.
#'
#' Two PODs per setting, 100,000 SNPs each, same simulation as
#' `_R/baypass-pod-inputs.R` except the Ω:
#'   base       source Ω (the current method at this size: the reference
#'              for FMD noise at 100,000 SNPs)
#'   corrected  Ω*
#' Each runs through the core model (`03-select-time/05-baypass-pod-test`),
#' and the test compares FMD(re-estimate, Ω) between the two. The
#' correction passes if FMD for `corrected` is below 0.5 and below `base`.
#'
#' Part 1 (inputs): writes `tany_dir("baypass-pod",
#' "baypass-pod-test-inputs.tar.gz")` (upload to
#' `/staging/l/lnell/dna/baypass/`). Part 2 (results; runs when all 6
#' outputs are in `tany_dir("baypass-pod")`). Both report to
#' `_data/baypass-pod-correct.txt`.
#'

library(tidyverse)
library(geigen)

#' Analyses the first 2026-11 rerun (18 pools; archived):
ARCHIVE_18POOLS <- TRUE
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

out_file <- "_data/baypass-pod-correct.txt"
sink(out_file, split = TRUE)

mp_tags <- c("mp050", "mp075", "mp100")
variants <- c("base", "corrected")
#' Globals used by `simulate_pod()`:
n_pod <- 100000L
n_batch <- 50000L

set.seed(1514460719)

fmd.dist <- function(mat1, mat2) {
    return(sqrt(sum(log(geigen(mat1, mat2)$values)**2)))
}
read_om <- function(tar, member) {
    tmp <- tempfile()
    untar(tar, files = member, exdir = tmp)
    as.matrix(read.table(file.path(tmp, member)))
}
#' Ω re-estimated from the full POD (core model):
om_re_full <- function(tag) {
    nm <- sprintf("baypass-pod-%s-core-run1", tag)
    read_om(tany_dir("baypass-pod", paste0(nm, ".tar")),
            file.path(nm, sprintf("tany-pod-%s-core-run1_mat_omega.out", tag)))
}



# =============================================================================*
# ---- part 1: inputs ----
# =============================================================================*


test_tar <- tany_dir("baypass-pod", "baypass-pod-test-inputs.tar.gz")

if (!file.exists(test_tar)) {
    cat("== Part 1: test PODs (", n_pod, " SNPs each)\n", sep = "")
    out_dir <- file.path(tempdir(), "baypass-pod-test-inputs")
    dir.create(out_dir)
    for (tag in mp_tags) {
        src <- get_pod_src(tag)
        om <- src$omega
        dd <- sqrt(diag(om) / diag(om_re_full(tag)))
        om_star <- om * outer(dd, dd)
        ev <- eigen(om_star, symmetric = TRUE, only.values = TRUE)$values
        cat(tag, ": Omega* diagonal / Omega ", paste(range(round(
            diag(om_star) / diag(om), 3)), collapse = "-"),
            "; smallest eigenvalue ", signif(min(ev), 3), "\n", sep = "")
        stopifnot(min(ev) > 0)
        for (v in variants) {
            s <- src
            if (v == "corrected") s$omega <- om_star
            pod <- simulate_pod(s)
            d <- file.path(out_dir, sprintf("baypass-podtest_%s_%s", tag, v))
            dir.create(d)
            g <- matrix(0L, n_pod, 2 * ncol(pod$ref))
            g[, seq(1, ncol(g), 2)] <- pod$ref
            g[, seq(2, ncol(g), 2)] <- pod$dep - pod$ref
            data.table::fwrite(as.data.frame(g), file.path(d, "tany-pod.genobaypass"),
                               sep = " ", col.names = FALSE)
            write_lines(paste(src$hap, collapse = " "),
                        file.path(d, "tany-pod.poolsize"))
            write.table(s$omega, file.path(d, "source_mat_omega.out"),
                        row.names = FALSE, col.names = FALSE)
        }
    }
    old <- setwd(out_dir)
    tar(test_tar, files = list.files(".", pattern = "^baypass-podtest_"),
        compression = "gzip", tar = "internal")
    setwd(old)
    cat("Wrote", test_tar, "\n")
}



# =============================================================================*
# ---- part 2: results ----
# =============================================================================*


res_name <- function(tag, v) sprintf("baypass-podtest-%s-%s", tag, v)
have <- file.exists(tany_dir("baypass-pod", paste0(
    outer(mp_tags, variants, res_name), ".tar")))
if (all(have)) {
    cat("\n== Part 2: Omega re-estimated from each test POD vs. the",
        "real-data Omega (subsample 1)\n")
    res <- crossing(tag = mp_tags, variant = variants) |>
        mutate(x = map2(tag, variant, \(t, v) {
            nm <- res_name(t, v)
            om_re <- read_om(tany_dir("baypass-pod", paste0(nm, ".tar")),
                             file.path(nm, sprintf("tany-podtest-%s-%s_mat_omega.out",
                                                   t, v)))
            om <- get_pod_src(t)$omega
            od <- upper.tri(om)
            r <- diag(om_re) / diag(om)
            tibble(fmd = fmd.dist(om_re, om),
                   diag_ratio_min = min(r), diag_ratio_max = max(r),
                   trace_ratio = sum(diag(om_re)) / sum(diag(om)),
                   offdiag_cor = cor(om_re[od], om[od]))
        })) |>
        unnest(x)
    res |>
        mutate(across(where(is.double), \(x) signif(x, 3))) |>
        as.data.frame() |>
        print()
    pass <- res |>
        select(tag, variant, fmd) |>
        pivot_wider(names_from = variant, values_from = fmd) |>
        mutate(passes = corrected < 0.5 & corrected < base)
    cat("\nCorrection passes (FMD < 0.5 and below base):\n")
    print(as.data.frame(pass))
    #' What FMD is made of: log generalized eigenvalues of (re-estimate,
    #' real Omega), sorted; FMD = sqrt(sum of their squares). Also the
    #' share of FMD^2 from the 3 largest |log| values, and FMD between the
    #' correlation matrices (scale removed).
    cat("\nLog generalized eigenvalues, re-estimate vs. real Omega:\n")
    crossing(tag = mp_tags, variant = variants) |>
        mutate(x = map2(tag, variant, \(t, v) {
            nm <- res_name(t, v)
            om_re <- read_om(tany_dir("baypass-pod", paste0(nm, ".tar")),
                             file.path(nm, sprintf("tany-podtest-%s-%s_mat_omega.out",
                                                   t, v)))
            om <- get_pod_src(t)$omega
            l <- sort(log(geigen(om_re, om)$values))
            tibble(min = min(l), q25 = quantile(l, 0.25), median = median(l),
                   q75 = quantile(l, 0.75), max = max(l),
                   top3_share = sum(sort(l^2, decreasing = TRUE)[1:3]) / sum(l^2),
                   fmd_cor = fmd.dist(cov2cor(om_re), cov2cor(om)))
        })) |>
        unnest(x) |>
        mutate(across(where(is.double), \(x) signif(x, 3))) |>
        as.data.frame() |>
        print()
} else {
    cat("\nPart 2: waiting for", sum(!have), "of", length(have),
        "test outputs\n")
}

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
