#'
#' Inputs for the ABC test (`_R/abc-design.md`): pools, trap counts, the
#' 2,000 observed 10-kb windows, read depth at every site of those windows
#' (used by both observed and simulated statistics), and observed
#' second-allele read counts.
#'
#' Windows: random among 10-kb windows with >= 5,000 unmasked sites in
#' all 18 pools. Unmasked counts come from the `length` column of the
#' `02-npstat` outputs; no other npstat column is read, so window choice is
#' blind to diversity.
#'
#' Alleles: at each site, the two bases with the most reads summed over
#' pools are the major and second allele; depth = major + second reads.
#' Reads of the other two bases (third and fourth) are dropped, and
#' counted to estimate the sequencing error rate: at a biallelic site,
#' errors go to each of the 3 other bases equally, so the share of reads
#' showing a third or fourth base is ~ (2/3) e.
#'
#' Reads `tany_dir("npstat")` and `tany_dir("snape")` (partially masked
#' gSYNC files; 18 SN pools). Writes `tany_dir("abc", "abc-inputs")` and
#' `tany_dir("abc", "abc-inputs.tar.gz")` (upload to
#' `/staging/l/lnell/dna/abc/`), and `_data/abc-inputs.txt` (checks).
#'

library(tidyverse)
library(data.table)

source("_R/data-dir.R")

out_txt <- "_data/abc-inputs.txt"
sink(out_txt, split = TRUE)

W <- 10000L
N_WIN <- 2000L
MIN_CALLABLE <- 5000L
MAX_DEP <- 255L

set.seed(1606617349)



# =============================================================================*
# ---- pools and trap counts ----
# =============================================================================*


pools <- read_lines("_bash/samples/time-pools.txt")
pool_df <- read_csv("_data/full-sample-info.csv", show_col_types = FALSE) |>
    distinct(biotech_id, n_adults) |>
    filter(biotech_id %in% pools) |>
    mutate(yy = as.integer(str_match(biotech_id, "^SN-[12]-(\\d{2})")[, 2]),
           year = ifelse(yy < 70, 2000L, 1900L) + yy,
           g12 = as.integer(str_match(biotech_id, "^SN-([12])-")[, 2]),
           gen = 2L * (year - 1977L) + g12 - 1L) |>
    arrange(gen) |>
    select(pool = biotech_id, gen, n_adults)
stopifnot(nrow(pool_df) == 18, !anyDuplicated(pool_df$gen),
          all(!is.na(pool_df$n_adults)))

trap <- read_csv("_data/tany-abundances.csv", show_col_types = FALSE) |>
    arrange(year, gen)
stopifnot(nrow(trap) == 78, trap$year[1] == 1977, trap$gen[1] == "a")

cat("== Pools (gen: 0 = 1977a, 77 = 2015b)\n")
print(as.data.frame(pool_df))



# =============================================================================*
# ---- windows ----
# =============================================================================*


#' Unmasked sites per 10-kb window (npstat `length`), all pools:
callable <- pool_df$pool |>
    map(\(p) {
        fread(cmd = paste0("gzip -dc '",
                           tany_dir("npstat", paste0(p, "_npstat_10kb.stat.gz")),
                           "'"),
              select = c("sequence", "window", "length")) |>
            as_tibble() |>
            mutate(pool = p)
    }) |>
    list_rbind()

ok_win <- callable |>
    group_by(sequence, window) |>
    summarize(n_pools = n(), min_len = min(length), .groups = "drop") |>
    filter(n_pools == nrow(pool_df), min_len >= MIN_CALLABLE)

cat("\nWindows with >= ", MIN_CALLABLE, " unmasked sites in all pools: ",
    nrow(ok_win), " of ", n_distinct(paste(callable$sequence,
                                          callable$window)), "\n", sep = "")

win_df <- ok_win |>
    slice_sample(n = N_WIN) |>
    mutate(contig = factor(sequence, levels = paste0("CONTIG", 1:45))) |>
    arrange(contig, window) |>
    mutate(locus = row_number(),
           start = (window - 1L) * W + 1L) |>
    select(locus, contig, window, start)

cat("Windows per contig:\n")
print(table(droplevels(win_df$contig)))



# =============================================================================*
# ---- reads at window sites ----
# =============================================================================*


#' Site counts (A, T, C, G) in the chosen windows for one pool, unmasked
#' sites only, as a data.table with `idx` = 0-based site index
#' ((locus - 1) * W + offset).
win_file <- tempfile(fileext = ".tsv")
win_df |>
    transmute(key = paste(contig, window, sep = ":"), locus) |>
    write_tsv(win_file, col_names = FALSE)

pool_counts <- function(p) {
    tar_file <- tany_dir("snape", paste0(p, "_snape.tar"))
    member <- sprintf("%s_snape/%s_snape_part_masked.sync.gz", p, p)
    awk <- paste0(
        "awk -v W=", W, " 'BEGIN{FS=OFS=\"\\t\"} ",
        "NR==FNR{loc[$1]=$2; next} ",
        "$4!=\".:.:.:.:.:.\"{w=int(($2-1)/W)+1; k=$1\":\"w; ",
        "if(k in loc){split($4,a,\":\"); ",
        "print (loc[k]-1)*W+($2-1)-(w-1)*W, a[1], a[2], a[3], a[4]}}' ",
        win_file, " -")
    cmd <- sprintf("tar -xOf '%s' '%s' | gzip -dc | %s", tar_file, member,
                   awk)
    fread(cmd = cmd, header = FALSE,
          col.names = c("idx", "A", "T", "C", "G"))
}

n_sites <- N_WIN * W
tot <- matrix(0L, n_sites, 4)
tmp_dir <- tempfile("abc-counts")
dir.create(tmp_dir)
for (j in seq_len(nrow(pool_df))) {
    x <- pool_counts(pool_df$pool[j])
    i <- x$idx + 1L
    tot[i, ] <- tot[i, ] + as.matrix(x[, c("A", "T", "C", "G")])
    saveRDS(x, file.path(tmp_dir, paste0(j, ".rds")), compress = FALSE)
    cat(pool_df$pool[j], ": ", nrow(x), " unmasked sites\n", sep = "")
}

#' Major and second base per site (by total reads over pools):
major <- max.col(tot, ties.method = "first")
tot[cbind(seq_len(n_sites), major)] <- -1L
second <- max.col(tot, ties.method = "first")
rm(tot)

depth <- array(as.raw(0), dim = c(W, nrow(pool_df), N_WIN))
alt_list <- vector("list", nrow(pool_df))
err <- tibble()
for (j in seq_len(nrow(pool_df))) {
    x <- readRDS(file.path(tmp_dir, paste0(j, ".rds")))
    i <- x$idx + 1L
    m <- as.matrix(x[, c("A", "T", "C", "G")])
    a1 <- m[cbind(seq_along(i), major[i])]
    a2 <- m[cbind(seq_along(i), second[i])]
    other <- rowSums(m) - a1 - a2
    d <- a1 + a2
    err <- bind_rows(err, tibble(pool = pool_df$pool[j],
                                 reads = sum(as.numeric(rowSums(m))),
                                 other = sum(as.numeric(other)),
                                 over_max = mean(d > MAX_DEP),
                                 median_depth = median(d)))
    keep <- d > 0L & d <= MAX_DEP
    locus <- x$idx[keep] %/% W + 1L
    offset <- x$idx[keep] %% W
    depth[cbind(offset + 1L, j, locus)] <- as.raw(d[keep])
    #' Second-allele counts where they can matter (>= 1 read here; sites
    #' are kept below only if some pool has >= 2):
    nz <- keep & a2 > 0L
    alt_list[[j]] <- data.table(idx = x$idx[nz], j = j, alt = a2[nz])
    rm(x, m)
}
unlink(tmp_dir, recursive = TRUE)

err <- err |>
    mutate(other_share = other / reads,
           e_est = 1.5 * other_share)
cat("\n== Third/fourth-base reads (sequencing error check; design e = 0.002)\n")
print(as.data.frame(err |> mutate(across(where(is.double), \(x) signif(x, 3)))))
cat("Pooled e estimate:", signif(1.5 * sum(err$other) / sum(err$reads), 3),
    "\n")

#' Sites where some pool has >= 2 second-allele reads (no other site can
#' be segregating, since a site segregates in a pool if its minor allele
#' has >= 2 reads):
alt_dt <- rbindlist(alt_list)
rm(alt_list)
idx2 <- sort(unique(alt_dt$idx[alt_dt$alt >= 2L]))
alt_dt <- alt_dt[idx %in% idx2]
alt_m <- matrix(0L, length(idx2), nrow(pool_df),
                dimnames = list(NULL, pool_df$pool))
alt_m[cbind(match(alt_dt$idx, idx2), alt_dt$j)] <- alt_dt$alt
obs <- data.table(locus = idx2 %/% W + 1L, offset = idx2 %% W, alt_m)
cat("\nObserved sites with >= 2 second-allele reads in some pool:",
    nrow(obs), "\n")



# =============================================================================*
# ---- write ----
# =============================================================================*


out_dir <- tany_dir("abc", "abc-inputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(pool_df, file.path(out_dir, "pools.tsv"))
write_lines(trap$N, file.path(out_dir, "trap.txt"))
write_tsv(win_df |> mutate(contig = as.character(contig)),
          file.path(out_dir, "windows.tsv"))
writeBin(as.vector(depth), file.path(out_dir, "depth.bin"))
fwrite(obs, file.path(out_dir, "obs-sites.tsv.gz"), sep = "\t")
file.copy(c("_bash/04-diversity/abc-priors.csv",
            "_bash/04-diversity/abc-priors-mu.csv"), out_dir, overwrite = TRUE)

tar_file <- tany_dir("abc", "abc-inputs.tar.gz")
old <- setwd(tany_dir("abc"))
tar(tar_file, files = "abc-inputs", compression = "gzip", tar = "internal")
setwd(old)
cat("\nWrote", out_dir, "and", tar_file, "\n")

sink()
read_lines(out_txt) |>
    str_trim("right") |>
    write_lines(out_txt)
