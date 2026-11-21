#'
#' Fine-scale diversity around the CONTIG18 locus (owner, 2026-11-20;
#' exploratory). Long-term balancing selection would raise diversity and
#' Tajima's D within a few kb of the selected site; npstat's 10-kb windows
#' can't resolve that (`_data/contig18-diversity.txt`).
#'
#' Per pool, read counts at every site of CONTIG18:323,000-423,000 (the
#' locus at 373,124-373,206 +/- ~50 kb) from the partially masked SNAPE
#' counts (`tany_dir("snape", "<pool>_snape.tar")`; masked sites dropped),
#' cached per pool in `_data/contig18-region/` (local). In 1-kb windows,
#' with the estimators of the ABC statistics (`_bash/bin/abc-stats.R`;
#' minimum minor-allele count 2, depth <= 255, the pool's two most common
#' bases): θ_π and θ_W per site, and θ_π / θ_W - 1 (sign of Tajima's D).
#' Pools: the 17 SN pools and the 18 spatial pools of 2019.
#'
#' Caveat: reads carrying the B allele map poorly, so the window holding
#' the locus (373 kb) reads low; flanking windows are the test.
#'
#' Writes `_data/contig18-finescale.txt` and `_plots/contig18-finescale.pdf`.
#'

library(tidyverse)

source("_R/data-dir.R")
source("_bash/bin/abc-stats.R")  # theta_denoms(), B_MIN, MAX_DEP

out_file <- "_data/contig18-finescale.txt"
sink(out_file, split = TRUE)

ctg <- "CONTIG18"
r_lo <- 323000
r_hi <- 423000
lo <- 373124
hi <- 373206
win <- 1000
cache_dir <- "_data/contig18-region"
dir.create(cache_dir, showWarnings = FALSE)

#' Pool sets to include (default both): `Rscript _R/contig18-finescale.R time`.
set_arg <- commandArgs(TRUE)
if (length(set_arg) == 0) set_arg <- c("time", "space")
sets <- list(time = read_lines("_bash/samples/time-pools.txt"),
             space = read_lines("_bash/samples/space-pools.txt"))[set_arg]
n_adults <- read_csv("_data/full-sample-info.csv", show_col_types = FALSE) |>
    distinct(biotech_id, n_adults) |>
    (\(x) setNames(x$n_adults, x$biotech_id))()



# =============================================================================*
# ---- extract (cached per pool) ----
# =============================================================================*


#' Sites of the region in one pool: position and counts of A, T, C, G.
#' The sync files list contigs in lexical order, so awk stops after the
#' region instead of reading the whole file.
region_counts <- function(p) {
    f <- file.path(cache_dir, paste0(p, ".rds"))
    if (file.exists(f)) return(read_rds(f))
    tar <- tany_dir("snape", paste0(p, "_snape.tar"))
    member <- sprintf("%s_snape/%s_snape_part_masked.sync.gz", p, p)
    awk <- sprintf(paste0("awk -v c=%s -v a=%d -v b=%d 'BEGIN{FS=OFS=\"\\t\"} ",
                          "$1==c && $2>=a && $2<=b && $4!=\".:.:.:.:.:.\"",
                          "{split($4,x,\":\"); print $2,x[1],x[2],x[3],x[4]} ",
                          "$1==c && $2>b{exit}'"), ctg, r_lo, r_hi)
    cmd <- sprintf("tar -xOf '%s' '%s' 2>/dev/null | gzip -dc 2>/dev/null | %s",
                   tar, member, awk)
    out <- suppressWarnings(system(cmd, intern = TRUE))
    x <- read_tsv(I(out), col_names = c("pos", "A", "T", "C", "G"),
                  col_types = "iiiii")
    stopifnot(nrow(x) > 1000)
    write_rds(x, f, compress = "gz")
    x
}



# =============================================================================*
# ---- statistics ----
# =============================================================================*


window_stats <- function(x, n) {
    m <- as.matrix(x[, c("A", "T", "C", "G")])
    s <- t(apply(m, 1, sort, decreasing = TRUE))
    c1 <- s[, 1] + s[, 2]
    m2 <- s[, 2]
    ok <- c1 >= 2 * B_MIN & c1 <= MAX_DEP
    dens <- theta_denoms(n)
    seg <- ok & m2 >= B_MIN
    h <- ifelse(seg, 2 * m2 * (c1 - m2) / (c1 * (c1 - 1)), 0)
    tibble(w = x$pos %/% win, ok = ok, seg = seg, h = h,
           dpi = ifelse(ok, dens$pi[pmin(pmax(c1, 1), MAX_DEP)], 0),
           dw = ifelse(ok, dens$w[pmin(pmax(c1, 1), MAX_DEP)], 0)) |>
        group_by(w) |>
        summarize(callable = sum(ok), S = sum(seg),
                  theta_pi = sum(h) / sum(dpi), theta_w = sum(seg) / sum(dw),
                  .groups = "drop") |>
        mutate(tdr = theta_pi / theta_w - 1)
}

res <- imap(sets, \(pools, set) {
    map(pools, \(p) {
        x <- region_counts(p)
        window_stats(x, 2L * n_adults[[p]]) |>
            mutate(pool = p, set = set)
    }) |>
        list_rbind()
}) |>
    list_rbind() |>
    filter(callable >= 500)

w_loc <- lo %/% win
summ <- res |>
    group_by(set, w) |>
    summarize(pools = n(), callable = median(callable),
              theta_pi = median(theta_pi), theta_w = median(theta_w),
              tdr = median(tdr, na.rm = TRUE), .groups = "drop") |>
    mutate(start_kb = w * win / 1000, dist_kb = (w - w_loc) * win / 1000)

cat("== Median over pools, 1-kb windows, ", ctg, ":", r_lo, "-", r_hi,
    " (locus window: ", w_loc * win / 1000, " kb)\n", sep = "")
for (st in names(sets)) {
    x <- summ |> filter(set == st)
    cat("\n-- ", st, " pools (", length(sets[[st]]), "): windows ",
        nrow(x), "; median over windows: theta_pi ",
        signif(median(x$theta_pi), 3), ", theta_pi / theta_W - 1 ",
        round(median(x$tdr), 3), "\n", sep = "")
    cat("Windows within 10 kb of the locus, with each window's rank among",
        "all windows of the region (1 = highest):\n")
    x |>
        mutate(pi_rank = rank(-theta_pi), tdr_rank = rank(-tdr)) |>
        filter(abs(dist_kb) <= 10) |>
        transmute(start_kb, dist_kb, callable = round(callable),
                  theta_pi = signif(theta_pi, 3), pi_rank,
                  tdr = round(tdr, 3), tdr_rank) |>
        as.data.frame() |>
        print()
    cat("Top 5 windows of the region by theta_pi:\n")
    x |>
        slice_max(theta_pi, n = 5) |>
        transmute(start_kb, dist_kb, theta_pi = signif(theta_pi, 3),
                  tdr = round(tdr, 3)) |>
        as.data.frame() |>
        print()
}

p1 <- summ |>
    pivot_longer(c(theta_pi, tdr), names_to = "stat") |>
    mutate(stat = recode(stat, theta_pi = "theta_pi (per site)",
                         tdr = "theta_pi / theta_W - 1")) |>
    ggplot(aes(start_kb + 0.5, value, color = set)) +
    geom_vline(xintercept = c(lo, hi) / 1000, linetype = 2, color = "gray50") +
    geom_line() +
    facet_wrap(~ stat, ncol = 1, scales = "free_y") +
    labs(x = paste(ctg, "position (kb)"), y = NULL, color = "Pools") +
    theme_bw()
pdf("_plots/contig18-finescale.pdf", width = 8, height = 6)
print(p1)
invisible(dev.off())

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
