#'
#' Do the pools with large excess divergence carry related individuals?
#' (Owner, 2026-11-19; see `_data/pool-excess.txt`.) A pool with few
#' independent genomes, e.g. many siblings, deviates from the population
#' along long shared chromosome stretches, so its deviations at linked SNPs
#' are correlated over long distances. Deviations from binomial sampling of
#' unrelated genomes are correlated only over the short distances of
#' population linkage disequilibrium.
#'
#' Per SN pool j and SNP i: deviation d_ij = p_ij - mean of the other SN
#' pools, standardized z_ij = d_ij / sqrt(e_i (1 - e_i) (c_ij + s_j)), with
#' e_i the SN mean frequency, c_ij the read and pool sampling variance
#' factor (1/depth + 1/n - 1/(depth n)), and s_j set so that mean z^2 = 1
#' in pool j. The sign of a deviation depends on which allele is counted,
#' so the statistic is squared: for SNP pairs (i, k),
#'   R = mean (z_i z_k)^2 - 1,
#' which is 0 for independent deviations and about 2 r^2 for deviations
#' correlated by r. R is computed by distance between the SNPs (pairs a
#' set number of SNPs apart on the same contig), and for unlinked pairs
#' (SNPs on different contigs, randomly paired) as the baseline. Then:
#' does R at long distances (10 kb to 1 Mb) rise with the pool's excess?
#'
#' Data: the cached counts of `_R/pool-admixture.R` (SN pools only).
#'
#' Writes `_data/pool-relatedness.txt` and `_plots/pool-relatedness.pdf`.
#'

library(tidyverse)

ARCHIVE_18POOLS <- TRUE
source("_R/data-dir.R")

out_file <- "_data/pool-relatedness.txt"
sink(out_file, split = TRUE)

set.seed(1145287661)
sn <- read_lines("_bash/samples/time-pools.txt")
dat <- read_rds(tany_dir("archive-18pools", "snape", "pool-admixture-counts.rds"))
hap <- read_csv("_data/full-sample-info.csv", show_col_types = FALSE) |>
    distinct(biotech_id, n_adults) |>
    filter(biotech_id %in% sn) |>
    (\(x) setNames(2 * x$n_adults, x$biotech_id))()
hap <- hap[sn]

keep <- rowSums(dat$dep[, sn] >= 10) == length(sn)
o <- order(factor(dat$contig[keep], levels = paste0("CONTIG", 1:45)),
           dat$pos[keep])
contig <- dat$contig[keep][o]
pos <- dat$pos[keep][o]
p <- (dat$ref[keep, sn] / dat$dep[keep, sn])[o, ]
dp <- dat$dep[keep, sn][o, ]
rm(dat)
e <- rowMeans(p)
v <- e * (1 - e)
ok <- v > 0.01
contig <- contig[ok]; pos <- pos[ok]; p <- p[ok, ]; dp <- dp[ok, ]; e <- e[ok]
v <- v[ok]
cat("SNPs:", nrow(p), "(depth >= 10 in all", length(sn), "SN pools; SN mean",
    "frequency 0.01-0.99)\n")

z <- sapply(seq_along(sn), \(j) {
    d <- p[, j] - rowMeans(p[, -j])
    cc <- 1 / dp[, j] + 1 / hap[j] - 1 / (dp[, j] * hap[j])
    #' s_j so that mean z^2 = 1:
    s <- uniroot(\(s) mean(d^2 / (v * (cc + s))) - 1, c(0, 10))$root
    d / sqrt(v * (cc + s))
})
colnames(z) <- sn

offsets <- c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000)
pairs <- map(offsets, \(off) {
    a <- seq_len(nrow(z) - off)
    a <- a[contig[a] == contig[a + off]]
    b <- a + off
    tibble(i = a, k = b, dist = pos[b] - pos[a])
}) |>
    list_rbind() |>
    mutate(bin = cut(dist, c(0, 100, 1e3, 1e4, 1e5, 1e6, Inf),
                     labels = c("<100 bp", "0.1-1 kb", "1-10 kb", "10-100 kb",
                                "0.1-1 Mb", ">1 Mb")))
#' Unlinked: random pairs on different contigs.
n_un <- 2e6
ui <- sample.int(nrow(z), n_un, replace = TRUE)
uk <- sample.int(nrow(z), n_un, replace = TRUE)
un <- contig[ui] != contig[uk]
ui <- ui[un]
uk <- uk[un]

r_tab <- map(seq_along(sn), \(j) {
    zz <- z[, j]
    linked <- pairs |>
        mutate(r = (zz[i] * zz[k])^2) |>
        group_by(bin) |>
        summarize(R = mean(r) - 1, n = n(), .groups = "drop")
    bind_rows(linked,
              tibble(bin = "unlinked", R = mean((zz[ui] * zz[uk])^2) - 1,
                     n = length(ui))) |>
        mutate(pool = sn[j])
}) |>
    list_rbind()

ex <- read_lines("_data/pool-excess.txt")
ex <- str_match(ex, "^\\s*\\d+\\s+(SN-\\S+)\\s+\\d+\\s+\\d+\\s+[\\d.]+\\s+([\\d.e-]+)")
ex <- tibble(pool = ex[, 2], excess = as.numeric(ex[, 3])) |> filter(!is.na(pool))

cat("\nPairs per distance bin:\n")
print(as.data.frame(r_tab |> filter(pool == sn[1]) |> select(bin, n)))
cat("\n== R = mean (z_i z_k)^2 - 1 by distance and pool (sorted by excess)\n")
r_tab |>
    left_join(ex, by = "pool") |>
    mutate(R = round(R, 3)) |>
    select(pool, excess, bin, R) |>
    pivot_wider(names_from = bin, values_from = R) |>
    arrange(desc(excess)) |>
    as.data.frame() |>
    print()

cat("\n== Spearman correlation of R (minus the unlinked baseline) with the",
    "per-pool excess, by distance\n")
r_tab |>
    group_by(pool) |>
    mutate(R_adj = R - R[bin == "unlinked"]) |>
    ungroup() |>
    filter(bin != "unlinked") |>
    left_join(ex, by = "pool") |>
    group_by(bin) |>
    summarize(rho = round(cor(R_adj, excess, method = "spearman"), 3),
              p = signif(cor.test(R_adj, excess, method = "spearman",
                                  exact = FALSE)$p.value, 2),
              .groups = "drop") |>
    as.data.frame() |>
    print()

p1 <- r_tab |>
    filter(bin != "unlinked") |>
    left_join(ex, by = "pool") |>
    ggplot(aes(bin, R, group = pool, color = excess)) +
    geom_line() +
    geom_point() +
    scale_color_viridis_c() +
    labs(x = "Distance between SNPs", y = "R = mean (z_i z_k)^2 - 1",
         color = "Excess") +
    theme_classic()
pdf("_plots/pool-relatedness.pdf", width = 8, height = 5)
print(p1)
invisible(dev.off())

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
