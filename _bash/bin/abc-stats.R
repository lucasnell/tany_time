#!/usr/bin/env Rscript

#'
#' Summary statistics for the ABC test (`_R/abc-design.md`), computed by
#' the same code for observed and simulated data.
#'
#' Usage:
#'   Rscript abc-stats.R observed <inputs dir> <out.tsv>
#'   Rscript abc-stats.R sims <inputs dir> <model> <first seed> <n sims> \
#'       <out.tsv> [extra arguments for abc-sim.py]
#'
#' Inputs dir (from `_R/abc-inputs.R`):
#'   pools.tsv          pool, gen (0 = 1977a ... 77 = 2015b), n_adults
#'   trap.txt           trap count per generation (78 lines)
#'   windows.tsv        locus (1-based), contig, window, start
#'   depth.bin          uint8 read depth (major + second allele) per site,
#'                      pool, and locus (site fastest, then pool, then
#'                      locus); 0 = masked or depth > 255
#'   obs-sites.tsv.gz   observed sites with >= 2 second-allele reads in
#'                      some pool: locus, offset (0-based), then one
#'                      column per pool with the second-allele read count
#'   abc-priors.csv     priors for `abc-sim.py`
#'
#' `sims` calls `abc-sim.py` (next to this script, or `ABC_SIM_PY`) with
#' the python in `ABC_PYTHON` (default `python3`), then adds pool-seq
#' sampling: reads at each site ~ Binomial(depth, p (1 - e) + (1 - p) e / 3),
#' with the observed depth at that site, pool, and locus, and e the
#' sequencing error rate per read (0.00023, estimated from third- and
#' fourth-base reads in `_data/abc-inputs.txt`; owner, 2026-11-12;
#' `ABC_ERR` overrides). Sites that don't segregate in the sample
#' but show >= 2 error reads of one base are added too (Poisson thinning
#' over callable sites), since the observed data have them.
#'
#' Statistics (258 columns; see `_R/abc-design.md`):
#'   pi_<pool>, thw_<pool>   θ estimates per site, unbiased under the
#'                           standard neutral model given each site's
#'                           depth, the pool's haploid size, and the
#'                           minimum minor-allele count b = 2 (the
#'                           correction principle of Kofler et al. 2011)
#'   tdr_<pool>              pi / thw - 1 (the sign of Tajima's D)
#'   rare_<pool>             share of segregating sites with minor read
#'                           frequency < 0.05
#'   fc_<pool1>_<pool2>      temporal F between consecutive pools (17) and
#'                           from the first pool to each later one except
#'                           the second (16): Jónás et al. 2016 sampling
#'                           plan II estimator (poolSeq `estimateNe`,
#'                           method "P.planII"), on sites segregating in
#'                           either pool
#'   fst_<pool1>_<pool2>     pairwise F_ST (153), poolfstat
#'                           `compute.pairwiseFST` (Anova; Hivert et al.
#'                           2018)
#'

suppressPackageStartupMessages({
    library(stats)
    library(utils)
})

W <- 10000L      # locus / window length
B_MIN <- 2L      # minimum minor-allele read count for a segregating site
RARE <- 0.05     # rare-allele threshold (minor read frequency)
MAX_DEP <- 255L  # depths above this are masked (stored as uint8)



# =============================================================================*
# ---- setup ----
# =============================================================================*


#' Read the inputs and precompute what doesn't change among simulations.
abc_setup <- function(in_dir) {
    pools <- read.delim(file.path(in_dir, "pools.tsv"))
    pools <- pools[order(pools$gen), ]
    windows <- read.delim(file.path(in_dir, "windows.tsv"))
    L <- nrow(windows)
    P <- nrow(pools)
    depth <- readBin(file.path(in_dir, "depth.bin"), what = "raw",
                     n = W * P * L)
    stopifnot(length(depth) == W * P * L)
    n_hap <- 2L * pools$n_adults
    #' Callable sites per pool by depth (1 to MAX_DEP):
    hist <- vapply(seq_len(P), \(j) {
        d <- as.integer(depth[pool_index(j, P, L)])
        tabulate(d, nbins = MAX_DEP)
    }, numeric(MAX_DEP))
    #' Denominators: expected π- and S-type sums per site under the
    #' standard neutral model, per unit θ.
    dens <- lapply(seq_len(P), \(j) theta_denoms(n_hap[j]))
    d_pi <- vapply(seq_len(P), \(j) sum(hist[, j] * dens[[j]]$pi), 0)
    d_w <- vapply(seq_len(P), \(j) sum(hist[, j] * dens[[j]]$w), 0)
    list(pools = pools, windows = windows, L = L, P = P, depth = depth,
         n_hap = n_hap, hist = hist, d_pi = d_pi, d_w = d_w,
         pairs_t = temporal_pairs(P), pairs_f = utils::combn(P, 2))
}

#' Indices into `depth` of every site of pool j.
pool_index <- function(j, P, L) {
    starts <- (j - 1L) * W + (seq_len(L) - 1L) * W * P
    rep(starts, each = W) + rep(seq_len(W), L)
}

#' For depths 1 to MAX_DEP and haploid pool size n: sum over derived
#' counts k (1 to n - 1) of (1 / k) x E[statistic | depth, k / n], with the
#' statistic restricted to b <= m <= c - b.
theta_denoms <- function(n) {
    k <- seq_len(n - 1L)
    out <- lapply(seq_len(MAX_DEP), \(cc) {
        if (cc < 2L * B_MIN) return(c(pi = 0, w = 0))
        m <- B_MIN:(cc - B_MIN)
        pm <- outer(m, k / n, \(m, p) dbinom(m, cc, p))  # m x k
        h <- 2 * m * (cc - m) / (cc * (cc - 1))
        c(pi = sum(colSums(h * pm) / k), w = sum(colSums(pm) / k))
    })
    out <- do.call(rbind, out)
    list(pi = out[, "pi"], w = out[, "w"])
}

#' Temporal pairs (pools in time order): consecutive pairs, then the first
#' pool with each later pool except the second (already consecutive).
temporal_pairs <- function(P) {
    cbind(rbind(1:(P - 1), 2:P), rbind(1L, 3:P))
}



# =============================================================================*
# ---- statistics ----
# =============================================================================*


#' alt, dep: sites x pools integer matrices (pools in time order; dep 0 =
#' missing). Sites need only include those where some pool has >= 2
#' second-allele reads; no other site adds to any statistic.
compute_stats <- function(alt, dep, setup) {
    P <- setup$P
    nm <- setup$pools$pool
    #' Per pool, one column at a time (memory: simulations can have > 10
    #' million sites):
    seg <- matrix(FALSE, nrow(alt), P)
    pi_ <- thw <- rare <- numeric(P)
    for (j in seq_len(P)) {
        d <- dep[, j]
        m <- pmin(alt[, j], d - alt[, j])
        sj <- d > 0L & m >= B_MIN
        seg[, j] <- sj
        pi_[j] <- sum(2 * m[sj] * (d[sj] - m[sj]) / (d[sj] * (d[sj] - 1))) /
            setup$d_pi[j]
        thw[j] <- sum(sj) / setup$d_w[j]
        rare[j] <- sum(m[sj] < RARE * d[sj]) / sum(sj)
    }
    rm(d, m, sj)
    out <- c(setNames(pi_, paste0("pi_", nm)),
             setNames(thw, paste0("thw_", nm)),
             setNames(pi_ / thw - 1, paste0("tdr_", nm)),
             setNames(rare, paste0("rare_", nm)))

    #' Temporal F (Jónás et al. 2016, plan II), on sites covered in both
    #' pools and segregating in either. Jónás et al. keep sites polymorphic
    #' in the first sample instead, which with no drift gave F of -0.012
    #' to -0.0007 here; the rule used gives +0.0038 at 15x/30x depth,
    #' +0.0021 at 30x, and +0.0001 at 100x (`_bash/tests/abc-stats-check.R`
    #' → `_data/abc-stats-check.txt`). The bias is the same for observed
    #' and simulated data, so it doesn't affect the ABC, but F converted
    #' to Ne at low depth is biased.
    fc <- apply(setup$pairs_t, 2, \(ab) {
        a <- ab[1]; b <- ab[2]
        ok <- dep[, a] > 0L & dep[, b] > 0L & (seg[, a] | seg[, b])
        x <- alt[ok, a] / dep[ok, a]
        y <- alt[ok, b] / dep[ok, b]
        z <- (x + y) / 2
        s0 <- dep[ok, a]
        st <- dep[ok, b]
        n0 <- setup$n_hap[a]
        nt <- setup$n_hap[b]
        c0 <- 1 / s0 + 1 / n0 - 1 / (s0 * n0)
        ct <- 1 / st + 1 / nt - 1 / (st * nt)
        sum((x - y)^2 - (z - x * y) * (c0 + ct)) /
            sum((z - x * y) * (1 - ct))
    })
    names(fc) <- paste("fc", nm[setup$pairs_t[1, ]], nm[setup$pairs_t[2, ]],
                       sep = "_")
    out <- c(out, fc)

    #' Pairwise F_ST (poolfstat):
    out <- c(out, pairwise_fst(alt, dep, seg, setup))
    out
}

pairwise_fst <- function(alt, dep, seg, setup) {
    nm <- setup$pools$pool
    pr <- setup$pairs_f
    lab <- paste("fst", nm[pr[1, ]], nm[pr[2, ]], sep = "_")
    if (!requireNamespace("poolfstat", quietly = TRUE)) {
        warning("poolfstat not installed: F_ST set to NA")
        return(setNames(rep(NA_real_, ncol(pr)), lab))
    }
    keep <- rowSums(seg) > 0L
    pd <- new("pooldata")
    pd@npools <- setup$P
    pd@nsnp <- sum(keep)
    pd@refallele.readcount <- alt[keep, , drop = FALSE]
    pd@readcoverage <- dep[keep, , drop = FALSE]
    pd@snp.info <- data.frame(Chromosome = "1", Position = seq_len(sum(keep)),
                              RefAllele = NA, AltAllele = NA)
    pd@poolsizes <- setup$n_hap
    pd@poolnames <- nm
    colnames(pd@refallele.readcount) <- colnames(pd@readcoverage) <- nm
    f <- poolfstat::compute.pairwiseFST(pd, min.cov.per.pool = 1L,
                                        verbose = FALSE)
    fm <- f@PairwiseFSTmatrix
    setNames(fm[cbind(pr[1, ], pr[2, ])], lab)
}



# =============================================================================*
# ---- observed ----
# =============================================================================*


observed_stats <- function(in_dir, setup) {
    obs <- read.delim(gzfile(file.path(in_dir, "obs-sites.tsv.gz")))
    alt <- as.matrix(obs[, setup$pools$pool])
    storage.mode(alt) <- "integer"
    dep <- site_depths(obs$locus, obs$offset, setup)
    alt[dep == 0L] <- 0L
    compute_stats(alt, dep, setup)
}

#' Depths (sites x pools) at loci (1-based) and offsets (0-based).
site_depths <- function(locus, offset, setup) {
    P <- setup$P
    base <- offset + (locus - 1L) * W * P + 1L
    vapply(seq_len(P), \(j) as.integer(setup$depth[base + (j - 1L) * W]),
           integer(length(base)))
}



# =============================================================================*
# ---- simulations ----
# =============================================================================*


#' P(the most common of 3 error bases has >= 2 reads) at depth c, error
#' rate e (errors equally likely to each other base).
p_err2 <- function(cc, e) {
    q <- e / 3
    #' P(each of the 3 error bases has <= 1 read): multinomial terms with
    #' j bases at 1 read and 3 - j at 0.
    p_ok <- vapply(cc, \(x) {
        j <- 0:min(3L, x)
        sum(choose(3, j) * exp(lfactorial(x) - lfactorial(x - j)) *
                q^j * (1 - e)^(x - j))
    }, 0)
    pmax(1 - p_ok, 0)
}

#' Error-only sites with >= 2 reads of one error base, per pool: Poisson
#' thinning over all sites (rate p_err2(depth) at callable sites).
error_sites <- function(setup, e) {
    P <- setup$P
    L <- setup$L
    q2 <- c(0, p_err2(seq_len(MAX_DEP), e))
    qmax <- max(q2)
    if (qmax == 0) {
        return(data.frame(locus = integer(0), offset = integer(0),
                          pool = integer(0), alt = integer(0)))
    }
    q <- e / 3
    out <- lapply(seq_len(P), \(j) {
        K <- rpois(1, W * L * qmax)
        u <- sample.int(W * L, K, replace = TRUE) - 1L
        locus <- u %/% W + 1L
        offset <- u %% W
        d <- as.integer(setup$depth[offset + (locus - 1L) * W * P +
                                        (j - 1L) * W + 1L])
        acc <- runif(K) < q2[d + 1L] / qmax
        d <- d[acc]
        #' Approximately the count given >= 2:
        a <- 2L + rbinom(length(d), pmax(d - 2L, 0L), q)
        data.frame(locus = locus[acc], offset = offset[acc],
                   pool = rep(j, length(d)), alt = pmin(a, d))
    })
    do.call(rbind, out)
}

#' Pool-seq reads for one simulation's output prefix.
sim_reads <- function(prefix, setup, e) {
    P <- setup$P
    pos <- readBin(paste0(prefix, ".pos.bin"), "integer", n = 1e9, size = 4)
    cnt <- readBin(paste0(prefix, ".counts.bin"), "integer",
                   n = length(pos) * P, size = 2)
    cnt <- matrix(cnt, ncol = P, byrow = TRUE)
    locus <- pos %/% W + 1L
    offset <- pos %% W
    rm(pos)
    dep <- site_depths(locus, offset, setup)
    #' Reads, one pool at a time (memory):
    alt <- matrix(0L, nrow(cnt), P)
    for (j in seq_len(P)) {
        p <- cnt[, j] / setup$n_hap[j]
        alt[, j] <- rbinom(nrow(cnt), dep[, j], p * (1 - e) + (1 - p) * e / 3)
    }
    rm(cnt, p)
    #' Keep only sites where some pool has >= 2 minor-allele reads (no
    #' other site adds to any statistic):
    keep <- logical(nrow(alt))
    for (j in seq_len(P)) {
        keep <- keep | pmin(alt[, j], dep[, j] - alt[, j]) >= B_MIN
    }
    alt <- alt[keep, , drop = FALSE]
    dep <- dep[keep, , drop = FALSE]
    key <- (locus * W + offset)[keep]
    rm(locus, offset, keep)
    #' Error-only sites, unless the site already segregates in the sample
    #' (their minor allele has >= 2 reads by construction):
    es <- error_sites(setup, e)
    es <- es[!(es$locus * W + es$offset) %in% key, ]
    if (nrow(es) > 0) {
        ekey <- es$locus * W + es$offset
        u <- unique(ekey)
        i <- match(ekey, u)
        e_alt <- matrix(0L, length(u), P)
        e_alt[cbind(i, es$pool)] <- es$alt
        e_dep <- site_depths(u %/% W, u %% W, setup)
        alt <- rbind(alt, e_alt)
        dep <- rbind(dep, e_dep)
    }
    list(alt = alt, dep = dep)
}

run_sims <- function(in_dir, setup, model, seed0, n_sims, out_file,
                     extra = character(0)) {
    py <- Sys.getenv("ABC_PYTHON", "python3")
    e <- as.numeric(Sys.getenv("ABC_ERR", "0.00023"))
    sim_py <- Sys.getenv("ABC_SIM_PY",
                         file.path(dirname(this_file()), "abc-sim.py"))
    tmp <- tempfile("abc")
    write_lines <- function(x) cat(paste(x, collapse = "\t"), "\n", sep = "",
                                   file = out_file, append = TRUE)
    header_done <- file.exists(out_file) && file.size(out_file) > 0
    for (i in seq_len(n_sims) - 1L) {
        seed <- seed0 + i
        set.seed(seed)
        t0 <- proc.time()[["elapsed"]]
        st <- system2(py, c(sim_py, "--model", model, "--seed", seed,
                            "--out", tmp,
                            "--priors", file.path(in_dir, "abc-priors.csv"),
                            "--trap", file.path(in_dir, "trap.txt"),
                            "--pools", file.path(in_dir, "pools.tsv"),
                            "--n-loci", setup$L, extra))
        if (st != 0) stop("abc-sim.py failed for seed ", seed)
        t_sim <- proc.time()[["elapsed"]] - t0
        pars <- read.delim(paste0(tmp, ".params.tsv"))
        rd <- sim_reads(tmp, setup, e)
        s <- compute_stats(rd$alt, rd$dep, setup)
        t_all <- proc.time()[["elapsed"]] - t0
        row <- c(model = model, seed = seed,
                 setNames(pars$value, pars$param)[setdiff(pars$param,
                                                          c("model", "seed"))],
                 t_sim = round(t_sim, 1), t_all = round(t_all, 1),
                 setNames(format(s, digits = 6), names(s)))
        if (!header_done) {
            write_lines(names(row))
            header_done <- TRUE
        }
        write_lines(row)
        file.remove(paste0(tmp, c(".params.tsv", ".pos.bin", ".counts.bin")))
    }
}

this_file <- function() {
    a <- commandArgs(FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f) == 0) "abc-stats.R" else f
}



# =============================================================================*
# ---- main ----
# =============================================================================*


if (sys.nframe() == 0L) {
    args <- commandArgs(TRUE)
    mode <- args[1]
    in_dir <- args[2]
    setup <- abc_setup(in_dir)
    if (mode == "observed") {
        s <- observed_stats(in_dir, setup)
        write.table(data.frame(stat = names(s), value = s),
                    args[3], sep = "\t", quote = FALSE, row.names = FALSE)
    } else if (mode == "sims") {
        run_sims(in_dir, setup, model = args[3], seed0 = as.integer(args[4]),
                 n_sims = as.integer(args[5]), out_file = args[6],
                 extra = args[-(1:6)])
    } else {
        stop("mode must be 'observed' or 'sims'")
    }
}
