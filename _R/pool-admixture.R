#'
#' Are SN pools pulled toward other lakes? (Owner, 2026-11-19.) Test of
#' whether the per-pool excess divergence (`_data/pool-excess.txt`) is
#' partly immigrant ancestry. Sources named by the owner:
#'   very closely related lakes:   Mik-19_S15, Skj-19_S22, Vik-19_S41
#'   about as related as Mývatn sites to each other:  Sva-19_S39
#'   Mývatn sites, 2019:           MyBR-19_S17, MyKS-19-B_S18, MySN-19_S19
#'
#' Statistic, per SN pool j and source L: f4(j, H1; L, H2) =
#' mean over SNPs of (p_j - p_H1)(p_L - p_H2), with H1 and H2 two random
#' halves of the other SN pools (16; SN-1-95 excluded), averaged over 20
#' random splits. Using disjoint halves means no pool's sampling noise
#' enters both factors, so the expectation is 0 unless j shares ancestry
#' with L beyond what the other SN pools do. If a share q of pool j's
#' genomes came from L, f4 ~ q x Var(p_L - p_SN), so q ~ f4 / f4(L, H1;
#' L, H2) (the denominator is the same statistic with L in place of j:
#' L's own divergence from SN, free of L's sampling noise because the two
#' factors share only p_L... which they do, so the denominator subtracts
#' L's sampling variance). SE by block jackknife (1-Mb blocks). Then the
#' correlation of each source's q with the per-pool excess.
#'
#' Data: `count_masked_all_noblanks.sync.gz` of the first rerun's combined
#' counts (sites with data in all 36 pools; copied from the staging
#' archive to `tany_dir("archive-18pools", "snape")`), SNPs with pooled
#' minor-allele frequency >= 0.05 over the pools used.
#'
#' Writes `_data/pool-admixture.txt`.
#'

library(tidyverse)

ARCHIVE_18POOLS <- TRUE
source("_R/data-dir.R")

out_file <- "_data/pool-admixture.txt"
sink(out_file, split = TRUE)

set.seed(1979318307)
sn <- read_lines("_bash/samples/time-pools.txt")
sources <- c(Mik = "Mik-19_S15", Skj = "Skj-19_S22", Vik = "Vik-19_S41",
             Sva = "Sva-19_S39", MyBR = "MyBR-19_S17", MyKS = "MyKS-19-B_S18",
             MySN = "MySN-19_S19")
use <- c(sn, sources)



# =============================================================================*
# ---- counts ----
# =============================================================================*


#' Major-allele read count and depth (major + second allele, by reads
#' summed over the pools used) per pool and SNP; cached on the drive.
cache <- tany_dir("archive-18pools", "snape", "pool-admixture-counts.rds")
if (!file.exists(cache)) {
    tar <- tany_dir("archive-18pools", "snape", "count_masked_all.tar")
    nm <- system(sprintf("tar -xOf '%s' count_masked_all/count_masked_all.names.gz | gzip -dc",
                         tar), intern = TRUE)
    cols <- match(use, nm) + 3
    stopifnot(!anyNA(cols))
    awk <- paste0(
        "awk -v cols='", paste(cols, collapse = ","), "' -v minf=0.05 '",
        "BEGIN{n=split(cols,C,\",\")} {",
        "for(b=1;b<=4;b++)T[b]=0;",
        "for(i=1;i<=n;i++){split($C[i],a,\":\");for(b=1;b<=4;b++){X[i,b]=a[b];T[b]+=a[b]}}",
        "m1=1;for(b=2;b<=4;b++)if(T[b]>T[m1])m1=b;",
        "m2=(m1==1)?2:1;for(b=1;b<=4;b++)if(b!=m1&&T[b]>T[m2])m2=b;",
        "if(T[m1]+T[m2]==0||T[m2]/(T[m1]+T[m2])<minf)next;",
        "printf \"%s\\t%s\",$1,$2;",
        "for(i=1;i<=n;i++)printf \"\\t%d\\t%d\",X[i,m1],X[i,m1]+X[i,m2];",
        "printf \"\\n\"}'")
    cmd <- sprintf("tar -xOf '%s' count_masked_all/count_masked_all_noblanks.sync.gz | gzip -dc | %s",
                   tar, awk)
    x <- data.table::fread(cmd = cmd, header = FALSE)
    ref <- as.matrix(x[, seq(3, ncol(x), 2), with = FALSE])
    dep <- as.matrix(x[, seq(4, ncol(x), 2), with = FALSE])
    colnames(ref) <- colnames(dep) <- use
    write_rds(list(contig = x[[1]], pos = x[[2]], ref = ref, dep = dep),
              cache, compress = "gz")
}
dat <- read_rds(cache)
keep <- rowSums(dat$dep >= 10) == length(use)
p <- dat$ref[keep, ] / dat$dep[keep, ]
dp <- dat$dep[keep, ]
block <- paste(dat$contig[keep], dat$pos[keep] %/% 1e6)
cat("SNPs (pooled minor frequency >= 0.05, depth >= 10 in all", length(use),
    "pools):", nrow(p), "; 1-Mb blocks:", n_distinct(block), "\n")



# =============================================================================*
# ---- f4 and admixture share ----
# =============================================================================*


#' Per-block sums, so the jackknife only re-sums.
bsum <- function(v) tapply(v, block, sum)
nb <- n_distinct(block)
jack <- function(num_b, den_b = NULL) {
    tot_n <- sum(num_b)
    if (is.null(den_b)) {
        est <- tot_n / length(block)
        loo <- (tot_n - num_b) / (length(block) - tabulate(factor(block)))
    } else {
        est <- tot_n / sum(den_b)
        loo <- (tot_n - num_b) / (sum(den_b) - den_b)
    }
    se <- sqrt((nb - 1) / nb * sum((loo - mean(loo))^2))
    list(est = unname(est), se = unname(se))
}
n_split <- 20
splits <- map(seq_len(n_split), \(s) sample(rep(1:2, length.out = length(sn) - 1)))

#' Sampling variance of a pool's frequency (reads), to remove from L's
#' own divergence in the denominator.
samp_var <- function(L) p[, L] * (1 - p[, L]) / pmax(dp[, L] - 1, 1)

res <- crossing(pool = sn, source = names(sources)) |>
    mutate(x = map2(pool, source, \(j, s) {
        L <- sources[[s]]
        others <- setdiff(sn, j)
        num <- 0
        den <- 0
        for (sp in splits) {
            h1 <- rowMeans(p[, others[sp == 1], drop = FALSE])
            h2 <- rowMeans(p[, others[sp == 2], drop = FALSE])
            num <- num + (p[, j] - h1) * (p[, L] - h2)
            den <- den + (p[, L] - h1) * (p[, L] - h2) - samp_var(L)
        }
        num_b <- bsum(num / n_split)
        den_b <- bsum(den / n_split)
        jf <- jack(num_b)
        jq <- jack(num_b, den_b)
        tibble(f4 = jf$est, f4_se = jf$se, z = jf$est / jf$se,
               q = jq$est, q_se = jq$se)
    })) |>
    unnest(x)

excess <- read_lines("_data/pool-excess.txt")
ex <- str_match(excess, "^\\s*\\d+\\s+(SN-\\S+)\\s+\\d+\\s+\\d+\\s+[\\d.]+\\s+([\\d.e-]+)")
ex <- tibble(pool = ex[, 2], excess = as.numeric(ex[, 3])) |> filter(!is.na(pool))
stopifnot(nrow(ex) == length(sn))

cat("\n== Admixture share q (f4 ratio) of each source in each SN pool\n")
cat("(z = f4 / SE; |z| > 3 is the usual threshold for a significant f4)\n")
res |>
    left_join(ex, by = "pool") |>
    mutate(cell = sprintf("%.3f (%.1f)", q, z)) |>
    select(pool, excess, source, cell) |>
    pivot_wider(names_from = source, values_from = cell) |>
    arrange(desc(excess)) |>
    as.data.frame() |>
    print()

cat("\n== Per source: pools with z > 3, and Spearman correlation of q with",
    "the per-pool excess\n")
res |>
    left_join(ex, by = "pool") |>
    group_by(source) |>
    summarize(pools_z_gt3 = sum(z > 3), pools_z_lt_neg3 = sum(z < -3),
              mean_q = round(mean(q), 3),
              rho_q_excess = round(cor(q, excess, method = "spearman"), 3),
              .groups = "drop") |>
    as.data.frame() |>
    print()

#' The f4 above also measures nearness in time: all sources are from 2019,
#' so SN pools from later years share more drift with them. To remove
#' that, each other lake is compared with the same site in 2019
#' (MySN-19): f4(j, H; L, MySN-19) = mean (p_j - p_H)(p_L - p_MySN), with
#' H all other SN pools. Positive: pool j is closer to lake L than to SN
#' itself in 2019, as immigrant ancestry from L (relative to SN) would
#' make it. Both L and MySN-19 are 2019 samples, so time cancels.
cat("\n== Other lakes vs. SN in 2019: f4(j, other SN pools; L, MySN-19)",
    "x 1e3 (z)\n")
ref_pool <- sources[["MySN"]]
rel <- crossing(pool = sn, source = c("Mik", "Skj", "Vik", "Sva", "MyBR", "MyKS")) |>
    mutate(x = map2(pool, source, \(j, s) {
        L <- sources[[s]]
        h <- rowMeans(p[, setdiff(sn, j), drop = FALSE])
        jf <- jack(bsum((p[, j] - h) * (p[, L] - p[, ref_pool])))
        tibble(f4 = jf$est, z = jf$est / jf$se)
    })) |>
    unnest(x) |>
    left_join(ex, by = "pool")
rel |>
    mutate(cell = sprintf("%.3f (%.1f)", 1e3 * f4, z)) |>
    select(pool, excess, source, cell) |>
    pivot_wider(names_from = source, values_from = cell) |>
    arrange(desc(excess)) |>
    as.data.frame() |>
    print()
cat("\nPer lake: pools with z > 3 or < -3, and Spearman correlation of f4",
    "with the excess:\n")
rel |>
    group_by(source) |>
    summarize(z_gt3 = sum(z > 3), z_lt_neg3 = sum(z < -3),
              rho_excess = round(cor(f4, excess, method = "spearman"), 3),
              .groups = "drop") |>
    as.data.frame() |>
    print()

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
