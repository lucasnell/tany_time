#'
#' Exploratory look at the CONTIG18 locus (owner, 2026-11-19). Not a
#' candidate: it doesn't pass the pre-set rule (q 0.21-0.24 for
#' `hm_logN_6` in mp050 and mp075; `_data/baypass-check-candidates.txt`),
#' and it was singled out after seeing the results, so nothing here is
#' evidence of association. The aim is to see whether the signal looks
#' like biology or like an artifact.
#'
#' Locus: 9 SNPs at CONTIG18:373,124-373,206 beyond the POD null for
#' `hm_logN_6` in mp050 and mp075, at the BF ceiling in all 10 runs.
#'
#'   1. BayPass statistics of all tested SNPs within 10 kb, every setting
#'      and covariate (`_data/baypass-qvalues.rds`); whether the 9 SNPs
#'      were tested in mp100.
#'   2. SNP density: tested SNPs per kb around the locus vs. the genome
#'      (high density in a short stretch can mean collapsed paralogs).
#'   3. Read counts per pool at the 9 SNPs (BayPass inputs, mp050):
#'      reference-allele frequency, depth relative to each pool's median
#'      depth over all its SNPs (collapsed paralogs show high depth), and
#'      missing pools.
#'   4. Association: Pearson r of frequency with the covariate across
#'      pools, and the same leaving out each pool in turn (is one pool
#'      driving it?).
#'   5. Genes within 10 kb (`~/_data/Tgraci_genes.gff3`) and their
#'      functional annotation (`~/_data/midgenomes/go-terms/`).
#'
#' Writes `_data/contig18-locus.txt` and `_plots/contig18-locus.pdf`.
#'

library(tidyverse)

#' Analyses the first 2026-11 rerun (18 pools; archived):
ARCHIVE_18POOLS <- TRUE
source("_R/data-dir.R")

out_file <- "_data/contig18-locus.txt"
sink(out_file, split = TRUE)

ctg <- "CONTIG18"
lo <- 373124
hi <- 373206
flank <- 10000
pools <- read_lines("_bash/samples/time-pools-18.txt")



# =============================================================================*
# ---- 1. BayPass statistics nearby ----
# =============================================================================*


q_df <- read_rds("_data/archive-18pools/baypass-qvalues.rds")
near <- q_df |>
    filter(contig == ctg, pos >= lo - flank, pos <= hi + flank)
cat("== 1. Tested SNPs within 10 kb of ", ctg, ":", lo, "-", hi, "\n", sep = "")
near |>
    group_by(tag, covariate) |>
    summarize(snps = n(), in_locus = sum(pos >= lo & pos <= hi),
              beyond_null = sum(beyond_null),
              locus_min_q = signif(suppressWarnings(min(q[pos >= lo & pos <= hi])), 3),
              flank_max_bf = round(max(bf_db[pos < lo | pos > hi]), 1),
              flank_max_z = round(max(z[pos < lo | pos > hi]), 2),
              .groups = "drop") |>
    as.data.frame() |>
    print()
cat("\nSNPs in the locus, hm_logN_6:\n")
near |>
    filter(pos >= lo, pos <= hi, covariate == "hm_logN_6") |>
    select(tag, sub, pos, bf_db, z, n_runs_ceiling, p_emp, q) |>
    mutate(across(c(bf_db, z), \(x) round(x, 2)),
           across(c(p_emp, q), \(x) signif(x, 3))) |>
    arrange(tag, pos) |>
    as.data.frame() |>
    print()
cat("\nTop 10 SNPs within 10 kb by z, other covariates (mp050):\n")
near |>
    filter(tag == "mp050", covariate != "hm_logN_6") |>
    group_by(covariate) |>
    slice_max(z, n = 5) |>
    select(covariate, pos, bf_db, z, q) |>
    mutate(across(c(bf_db, z), \(x) round(x, 2)), q = signif(q, 3)) |>
    as.data.frame() |>
    print()



# =============================================================================*
# ---- 2. SNP density ----
# =============================================================================*


snps50 <- q_df |>
    filter(tag == "mp050", covariate == "hm_logN_6") |>
    select(contig, pos)
genome_len <- 91518485  # assembly length (legacy diversity scripts)
dens_genome <- nrow(snps50) / genome_len * 1000
win <- snps50 |>
    filter(contig == ctg) |>
    mutate(kb = floor(pos / 1000)) |>
    count(kb)
cat("\n== 2. SNPs tested in mp050 per kb: genome mean ",
    round(dens_genome, 2), "; CONTIG18 1-kb windows: median ",
    median(win$n), ", 99th percentile ", quantile(win$n, 0.99),
    "; window at 373 kb: ", win$n[win$kb == 373], "\n", sep = "")
cat("1-kb windows 368-378 kb:\n")
print(as.data.frame(win |> filter(kb >= 368, kb <= 378)))



# =============================================================================*
# ---- 3-4. read counts and association ----
# =============================================================================*


#' Counts at the locus SNPs and each pool's median depth over its SNPs,
#' from the mp050 BayPass inputs (subsamples 1-3).
inp <- tany_dir("baypass-inputs", "baypass-inputs_mp050_3.tar.gz")
tmp <- tempfile()
untar(inp, exdir = tmp)
d <- file.path(tmp, "baypass-inputs_mp050")
cnt <- map(1:3, \(k) {
    det <- read_table(file.path(d, paste0("tany.snpdet.sub", k)),
                      col_names = c("contig", "pos", "ref", "alt"),
                      col_types = "cdcc")
    g <- data.table::fread(file.path(d, paste0("tany.genobaypass.sub", k)),
                           header = FALSE) |>
        as.matrix()
    ref <- g[, seq(1, ncol(g), 2)]
    dep <- ref + g[, seq(2, ncol(g), 2)]
    dep_na <- dep
    dep_na[dep_na == 0] <- NA
    med <- apply(dep_na, 2, median, na.rm = TRUE)
    i <- which(det$contig == ctg & det$pos >= lo & det$pos <= hi)
    list(med = med,
         loc = tibble(pos = rep(det$pos[i], each = ncol(ref)),
                      pool = rep(pools, length(i)),
                      ref = as.vector(t(ref[i, , drop = FALSE])),
                      depth = as.vector(t(dep[i, , drop = FALSE]))))
})
unlink(tmp, recursive = TRUE)
pool_med <- colMeans(do.call(rbind, map(cnt, "med")))
loc <- map(cnt, "loc") |>
    list_rbind() |>
    mutate(freq = ifelse(depth > 0, ref / depth, NA_real_),
           rel_depth = depth / pool_med[match(pool, pools)])

cov_df <- tibble(pool = pools) |>
    mutate(hm_logN_6 = scan("_data/archive-18pools/covariates/hm_logN_6.txt",
                            quiet = TRUE),
           logN = scan("_data/archive-18pools/covariates/logN.txt", quiet = TRUE),
           yy = as.integer(str_match(pool, "^SN-[12]-(\\d{2})")[, 2]),
           year = ifelse(yy < 70, 2000L, 1900L) + yy,
           gen = as.integer(str_match(pool, "^SN-([12])-")[, 2]),
           time = year + (gen - 1) / 2)
loc <- loc |> left_join(cov_df, by = "pool")

cat("\n== 3. Read counts at the", n_distinct(loc$pos), "locus SNPs (mp050)\n")
cat("Depth relative to each pool's median depth (1 = typical):\n")
loc |>
    group_by(pool) |>
    summarize(median_depth = pool_med[match(pool[[1]], pools)],
              locus_depth = median(depth), rel_depth = round(median(rel_depth), 2),
              missing = sum(depth == 0), .groups = "drop") |>
    arrange(match(pool, pools)) |>
    as.data.frame() |>
    print()

cat("\nReference-allele frequency by pool (rows) and SNP (columns):\n")
loc |>
    select(pool, pos, freq) |>
    mutate(freq = round(freq, 2)) |>
    pivot_wider(names_from = pos, values_from = freq) |>
    left_join(cov_df |> select(pool, hm_logN_6), by = "pool") |>
    arrange(hm_logN_6) |>
    mutate(hm_logN_6 = round(hm_logN_6, 2)) |>
    as.data.frame() |>
    print()

#' Mean frequency over the locus SNPs (they move together: one haplotype),
#' after orienting each SNP so its frequencies correlate positively with
#' those of the first SNP (reference alleles can sit on either haplotype).
first <- loc |> filter(pos == min(pos)) |> select(pool, f1 = freq)
flip <- loc |>
    left_join(first, by = "pool") |>
    group_by(pos) |>
    summarize(flip = cor(freq, f1, use = "complete.obs") < 0, .groups = "drop")
loc <- loc |>
    left_join(flip, by = "pos") |>
    mutate(freq_o = ifelse(flip, 1 - freq, freq))
cat("\nSNPs flipped to orient them with", min(loc$pos), ":",
    paste(flip$pos[flip$flip], collapse = ", "), "\n")
hap <- loc |>
    group_by(pool, hm_logN_6, time) |>
    summarize(freq = mean(freq_o, na.rm = TRUE), .groups = "drop")
r_all <- cor(hap$freq, hap$hm_logN_6, use = "complete.obs")
r_loo <- map_dbl(pools, \(p) {
    x <- hap |> filter(pool != p)
    cor(x$freq, x$hm_logN_6, use = "complete.obs")
})
cat("\n== 4. Mean locus frequency (oriented; the allele of SNP", min(loc$pos),
    ") vs. hm_logN_6 across pools: Pearson r =",
    round(r_all, 3), "\n")
cat("Leaving out each pool:\n")
print(as.data.frame(tibble(left_out = pools, r = round(r_loo, 3)) |>
                        arrange(abs(r))))
#' Is it abundance or time? hm_logN_6 is high before the 2003-2005 crash
#' and low after it, so a shift between eras would also correlate with it.
hap <- hap |> mutate(after_2003 = time > 2003)
cat("Pearson r with year:", round(cor(hap$freq, hap$time, use = "complete.obs"), 3),
    "; mean frequency before 2003:", round(mean(hap$freq[!hap$after_2003], na.rm = TRUE), 3),
    ", after:", round(mean(hap$freq[hap$after_2003], na.rm = TRUE), 3), "\n")
r_within <- hap |>
    group_by(after_2003) |>
    summarize(n = sum(!is.na(freq)),
              r_hm = cor(freq, hm_logN_6, use = "complete.obs"), .groups = "drop")
cat("Pearson r with hm_logN_6 within each era:\n")
print(as.data.frame(r_within |> mutate(r_hm = round(r_hm, 3))))
cat("Pearson r of mean locus frequency with relative depth (reference bias?):",
    loc |>
        group_by(pool) |>
        summarize(freq = mean(freq_o, na.rm = TRUE), rd = median(rel_depth)) |>
        filter(!is.na(freq)) |>
        (\(x) round(cor(x$freq, x$rd), 3))(), "\n")
cat("\nMean locus frequency by pool, in time order:\n")
print(as.data.frame(hap |> arrange(time) |>
                        mutate(freq = round(freq, 2), hm_logN_6 = round(hm_logN_6, 2))))
cat("Correlation among the locus SNPs' frequencies across pools (min, median):",
    loc |>
        select(pool, pos, freq) |>
        pivot_wider(names_from = pos, values_from = freq) |>
        select(-pool) |>
        cor(use = "pairwise.complete.obs") |>
        (\(m) round(c(min(m), median(m[upper.tri(m)])), 3))(), "\n")



# =============================================================================*
# ---- 5. genes ----
# =============================================================================*


genes <- read_tsv("~/_data/Tgraci_genes.gff3", comment = "#",
                  col_names = c("contig", "source", "type", "start", "end",
                                "score", "strand", "phase", "attributes"),
                  col_types = "ccciicccc") |>
    filter(contig == ctg, end >= lo - flank, start <= hi + flank) |>
    mutate(transcript = str_match(attributes, "ID=([^;]+)")[, 2],
           contains_locus = start <= lo & end >= hi)
desc <- read_lines("~/_data/midgenomes/go-terms/Tgraci_go-terms.tsv")[-1] |>
    str_split("\t") |>
    map(\(x) tibble(transcript = x[[1]],
                    desc = paste(str_remove(x[startsWith(x, "description:")],
                                            "^description:"), collapse = "; "))) |>
    list_rbind()
cat("\n== 5. Annotated transcripts within 10 kb\n")
genes |>
    left_join(desc, by = "transcript") |>
    select(transcript, start, end, strand, contains_locus, desc) |>
    mutate(desc = str_trunc(desc, 90)) |>
    as.data.frame() |>
    print()



# =============================================================================*
# ---- plot ----
# =============================================================================*


p1 <- loc |>
    ggplot(aes(hm_logN_6, freq, group = factor(pos))) +
    geom_line(data = loc |> arrange(hm_logN_6), alpha = 0.3) +
    geom_point(aes(size = depth), alpha = 0.5) +
    labs(x = "hm_logN_6 (log harmonic mean N, 6 generations)",
         y = "Reference-allele frequency", size = "Depth") +
    theme_classic()
p2 <- hap |>
    ggplot(aes(time, freq)) +
    geom_line(color = "gray60") +
    geom_point(aes(color = hm_logN_6), size = 2) +
    scale_color_viridis_c() +
    labs(x = "Year", y = "Mean locus frequency", color = "hm_logN_6") +
    theme_classic()
dir.create("_plots", showWarnings = FALSE)
pdf("_plots/contig18-locus.pdf", width = 7, height = 7)
print(p1)
print(p2)
invisible(dev.off())

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
