#'
#' What is the B haplotype at the CONTIG18 locus, from raw reads? (Owner,
#' 2026-11-19; exploratory, see `_R/contig18-locus.R`.) The 9 "SNPs" in
#' 82 bp of an 84% AT region could be a small indel misaligned as
#' substitutions. Reads are anchored on a unique, non-variable 20-mer just
#' right of the locus (CONTIG18:373,247-373,266; the left side is
#' repetitive AT sequence) and the 90 bp before the anchor in each read
#' (373,157-373,246 if no indel; 7 of the 9 SNPs) are taken, on the reference strand.
#' Distinct sequences are counted per pool and compared with the assembly:
#' length difference (indel) and mismatches.
#'
#' Pools: all 17 time pools. Reads: `/Volumes/lan_images4/ill/dna/trimmed/`.
#'
#' Writes `_data/contig18-flank.txt` and `_data/contig18-flank-counts.csv`.
#'

library(tidyverse)
library(seqinr)

out_file <- "_data/contig18-flank.txt"
sink(out_file, split = TRUE)

ctg <- "CONTIG18"
anc_start <- 373247
anc_len <- 20L
up_len <- 90L
trim_dir <- "/Volumes/lan_images4/ill/dna/trimmed"
#' All 17 time pools (2026-11-19; the first run had the 6 below):
#' SN-2-97, SN-2-02, SN-1-97 (mapped B common) and SN-1-99, SN-2-15,
#' SN-2-08 (mapped B rare).
pools <- read_lines("_bash/samples/time-pools.txt")

asm <- read.fasta("~/_data/midgenomes/assemblies/Tgraci_assembly.fasta.gz",
                  seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
revcomp <- \(x) chartr("ACGT", "TGCA", paste(rev(strsplit(x, "")[[1]]),
                                             collapse = ""))
s <- toupper(asm[[ctg]][[1]])
anchor <- substr(s, anc_start, anc_start + anc_len - 1)
ref_up <- substr(s, anc_start - up_len, anc_start - 1)
#' Anchor occurrences in the whole assembly, both strands:
n_occ <- sum(map_int(asm, \(x) {
    x <- toupper(x[[1]])
    sum(gregexpr(anchor, x, fixed = TRUE)[[1]] > 0) +
        sum(gregexpr(revcomp(anchor), x, fixed = TRUE)[[1]] > 0)
}))
rm(asm)
cat("Anchor ", ctg, ":", anc_start, "-", anc_start + anc_len - 1, " ",
    anchor, "; occurrences in the assembly (both strands): ", n_occ, "\n",
    sep = "")
cat("Assembly, ", up_len, " bp before the anchor:\n", ref_up, "\n", sep = "")
stopifnot(n_occ == 1)

#' Reads with the anchor (or its reverse complement), oriented to the
#' reference strand; the up_len bp before the anchor (reads that don't reach
#' that far are dropped).
pat <- tempfile(fileext = ".txt")
write_lines(c(anchor, revcomp(anchor)), pat)
#' Anchored reads are cached (local only, `.rds` is git-ignored) so the
#' window can change without scanning the reads again.
cache <- "_data/contig18-flank-reads.rds"
raw <- if (file.exists(cache)) read_rds(cache) else tibble(pool = character(),
                                                           read = character())
todo <- setdiff(pools, raw$pool)
#' One pool at a time, saving the cache after each (so an interrupted run
#' resumes where it stopped).
for (p in todo) {
    tar <- file.path(trim_dir, paste0("trimmed_", p, ".tar"))
    members <- sprintf("trimmed_%s_L002_R%i_001.fastq.gz", p, 1:2)
    cmd <- sprintf("tar -xOf '%s' %s | gzip -dc | awk 'NR %% 4 == 2' | grep -F -f '%s'",
                   tar, paste(members, collapse = " "), pat)
    r <- system(cmd, intern = TRUE)
    stopifnot(is.null(attr(r, "status")))
    raw <- bind_rows(raw, tibble(pool = p,
                                 read = ifelse(grepl(anchor, r, fixed = TRUE), r,
                                               map_chr(r, revcomp))))
    write_rds(raw, cache, compress = "gz")
}
reads <- map(pools, \(p) {
    r <- raw$read[raw$pool == p]
    pos <- regexpr(anchor, r, fixed = TRUE)
    up <- ifelse(pos > up_len, substr(r, pos - up_len, pos - 1), NA_character_)
    tibble(pool = rep(p, sum(!is.na(up))), up = up[!is.na(up)])
}) |>
    list_rbind()
n_anchor <- raw |> dplyr::count(pool, name = "n_anchor")

#' Alignment of each distinct sequence to the assembly's (global, to see
#' indels): length difference and edit distance.
seqs <- reads |>
    dplyr::count(pool, up) |>
    group_by(up) |>
    mutate(total = sum(n)) |>
    ungroup()
top <- seqs |>
    distinct(up, total) |>
    arrange(desc(total)) |>
    slice_head(n = 8) |>
    mutate(id = paste0("s", row_number()),
           edit = map_int(up, \(u) as.integer(adist(u, ref_up))),
           is_assembly = up == ref_up)
cat("\n== Most common sequences (all pools),", up_len, "bp before the anchor\n")
for (i in seq_len(nrow(top))) {
    cat(top$id[i], " (", top$total[i], " reads; edit distance to assembly ",
        top$edit[i], if (top$is_assembly[i]) "; = assembly" else "", "):\n",
        top$up[i], "\n", sep = "")
}

#' Where the main non-assembly sequence differs: alignment with indels.
alt <- top |> filter(!is_assembly) |> slice_head(n = 1)
if (nrow(alt) == 1) {
    al <- attr(adist(alt$up, ref_up, counts = TRUE), "counts")
    cat("\nMost common non-assembly sequence vs. assembly: insertions ",
        al[, , "ins"], ", deletions ", al[, , "del"], ", substitutions ",
        al[, , "sub"], "\n", sep = "")
    tr <- attr(adist(alt$up, ref_up, counts = TRUE), "trafos")
    cat("Edit path (M match, S substitution, I insertion, D deletion):\n",
        tr, "\n", sep = "")
}

cat("\n== Reads per pool by sequence (s1-s8 as above; other = rest)\n")
reads |>
    left_join(top |> select(up, id), by = "up") |>
    mutate(id = replace_na(id, "other")) |>
    dplyr::count(pool, id) |>
    pivot_wider(names_from = id, values_from = n, values_fill = 0) |>
    right_join(n_anchor, by = "pool") |>
    arrange(match(pool, pools)) |>
    as.data.frame() |>
    print()

#' Haplotype of every read: closer (edit distance) to the assembly (A) or
#' to the most common non-assembly sequence (B, s2 above); reads more than
#' 10 edits from both are "other". Raw B frequency per pool (compared with
#' the mapped frequency in `_R/contig18-raw-assoc.R`).
b_seq <- top$up[!top$is_assembly][1]
cls <- reads |>
    mutate(dA = as.vector(adist(up, ref_up)), dB = as.vector(adist(up, b_seq)),
           hap = case_when(pmin(dA, dB) > 10 ~ "other",
                           dA <= dB ~ "A", TRUE ~ "B")) |>
    dplyr::count(pool, hap) |>
    pivot_wider(names_from = hap, values_from = n, values_fill = 0) |>
    mutate(raw_B = round(B / (A + B), 2))
cat("\n== Reads by haplotype (A = assembly, B = s2; other > 10 edits from",
    "both)\n")
print(as.data.frame(cls |> arrange(match(pool, pools))))
#' For `_R/contig18-raw-assoc.R`:
write_csv(cls, "_data/contig18-flank-counts.csv")

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
