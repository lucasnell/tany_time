#'
#' Check of `_bash/bin/abc-stats.R` (ABC summary statistics) on pools with
#' known diversity, run locally from the repo root:
#'   Rscript _bash/tests/abc-stats-check.R
#' Output: `_data/abc-stats-check.txt`.
#'
#' Sites: population derived-allele frequency p with density ∝ θ / p on
#' (1e-6, 1) (standard neutral spectrum), so each pool's sample spectrum is
#' ~ θ / k. All pools sample the same population (no drift between them),
#' so temporal F should be ~ 0. Pool sizes (haploid) and depths span the
#' real ranges (28-80; 15-200x). With no sequencing error, π and θ_W
#' estimates should be ~ θ, and pi / thw - 1 ~ 0. Temporal F with no drift
#' is also compared between the site rule used (segregating in either
#' pool) and Jónás et al.'s (polymorphic in the first sample).
#' F_ST needs poolfstat (not installed locally) and is NA here.
#'

source("_bash/bin/abc-stats.R")

out_file <- "_data/abc-stats-check.txt"
sink(out_file, split = TRUE)

set.seed(1940229367)

theta <- 0.01
P <- 6L
L <- 200L
in_dir <- tempfile("abc-check")
dir.create(in_dir)
pools <- data.frame(pool = paste0("p", 1:P), gen = c(30, 35, 44, 50, 61, 77),
                    n_adults = c(14, 19, 24, 40, 40, 40))
write.table(pools, file.path(in_dir, "pools.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(locus = 1:L, contig = "C1", window = 1:L,
                       start = (0:(L - 1)) * W + 1),
            file.path(in_dir, "windows.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
#' Depths: pool mean 15 to 200, Poisson around it; ~20% masked.
mean_dep <- c(15, 30, 60, 100, 150, 200)
depth <- array(0L, c(W, P, L))
for (j in 1:P) {
    d <- rpois(W * L, mean_dep[j])
    d[runif(W * L) < 0.2 | d < 10] <- 0L
    depth[, j, ] <- pmin(d, MAX_DEP)
}
writeBin(as.raw(as.vector(depth)), file.path(in_dir, "depth.bin"))

setup <- abc_setup(in_dir)

#' Sites: positions with a population polymorphism (rate θ ln(1e6)).
n_poly <- rpois(1, theta * log(1e6) * W * L)
u <- sample.int(W * L, n_poly) - 1L
p <- exp(log(1e-6) * (1 - runif(n_poly)))
cnt <- sapply(setup$n_hap, \(n) rbinom(n_poly, n, p))

reads <- function(cnt, u, e) {
    dep <- site_depths(u %/% W + 1L, u %% W, setup)
    q <- sweep(cnt, 2, setup$n_hap, "/")
    q <- q * (1 - e) + (1 - q) * e / 3
    alt <- matrix(rbinom(length(dep), dep, q), ncol = P)
    list(alt = alt, dep = dep)
}

cat("theta =", theta, "; pools:", paste(setup$n_hap, collapse = ", "),
    "haploid; mean depths:", paste(mean_dep, collapse = ", "), "\n")

rd <- reads(cnt, u, e = 0)
s <- compute_stats(rd$alt, rd$dep, setup)
cat("\n== No sequencing error\n")
print(round(s[grepl("^(pi|thw|tdr|rare)_", names(s))], 4))
print(round(s[grepl("^fc_", names(s))], 4))

#' With errors (e = 0.002), including error-only sites as in `sim_reads`:
e <- 0.002
rd <- reads(cnt, u, e = e)
es <- error_sites(setup, e)
ekey <- es$locus * W + es$offset
key <- (u %/% W + 1L) * W + u %% W
es <- es[!ekey %in% key, ]
uu <- unique(es$locus * W + es$offset)
e_alt <- matrix(0L, length(uu), P)
e_alt[cbind(match(es$locus * W + es$offset, uu), es$pool)] <- es$alt
alt <- rbind(rd$alt, e_alt)
dep <- rbind(rd$dep, site_depths(uu %/% W, uu %% W, setup))
s2 <- compute_stats(alt, dep, setup)
cat("\n== Sequencing error e = 0.002 (bias expected; same in observed data)\n")
print(round(s2[grepl("^(pi|thw|tdr|rare)_", names(s2))], 4))
cat("Error-only sites added:", nrow(es), "\n")

#' Temporal F with no drift, two pools, 8 replicates per setting: the
#' site rule used (segregating in either pool) vs. Jónás et al.'s
#' (polymorphic in the first sample).
fc2 <- function(alt, dep, n_hap, ok) {
    x <- alt[ok, 1] / dep[ok, 1]; y <- alt[ok, 2] / dep[ok, 2]; z <- (x + y) / 2
    s0 <- dep[ok, 1]; st <- dep[ok, 2]
    c0 <- 1 / s0 + 1 / n_hap[1] - 1 / (s0 * n_hap[1])
    ct <- 1 / st + 1 / n_hap[2] - 1 / (st * n_hap[2])
    sum((x - y)^2 - (z - x * y) * (c0 + ct)) / sum((z - x * y) * (1 - ct))
}
fc_rep <- function(n_hap, md, seed) {
    set.seed(seed)
    n <- rpois(1, theta * log(1e6) * 2e6)
    p <- exp(log(1e-6) * (1 - runif(n)))
    cnt <- sapply(n_hap, \(nn) rbinom(n, nn, p))
    dep <- sapply(md, \(m) pmax(rpois(n, m), 10L))
    alt <- matrix(rbinom(length(dep), dep, sweep(cnt, 2, n_hap, "/")), ncol = 2)
    seg <- pmin(alt, dep - alt) >= B_MIN
    c(either = fc2(alt, dep, n_hap, seg[, 1] | seg[, 2]),
      first_poly = fc2(alt, dep, n_hap, alt[, 1] > 0 & alt[, 1] < dep[, 1]))
}
cat("\n== Temporal F with no drift (mean and SE of 8 replicates)\n")
for (cfg in list(list(c(28, 38), c(15, 30)), list(c(80, 80), c(30, 30)),
                 list(c(80, 80), c(100, 100)), list(c(28, 80), c(15, 100)))) {
    v <- sapply(1:8, \(s) fc_rep(cfg[[1]], cfg[[2]], s))
    cat("haploid ", paste(cfg[[1]], collapse = "/"), ", depth ",
        paste(cfg[[2]], collapse = "/"), ": segregating in either ",
        signif(mean(v[1, ]), 2), " (SE ", signif(sd(v[1, ]) / sqrt(8), 2),
        "); polymorphic in first ", signif(mean(v[2, ]), 2), " (SE ",
        signif(sd(v[2, ]) / sqrt(8), 2), ")\n", sep = "")
}

sink()
#' Printed vectors pad lines with spaces; strip them.
writeLines(sub("[[:space:]]+$", "", readLines(out_file)), out_file)
