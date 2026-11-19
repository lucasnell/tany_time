#'
#' What is the B haplotype at the CONTIG18 locus, from raw reads? (Owner,
#' 2026-11-19; exploratory, see `_R/contig18-locus.R`.) The 9 "SNPs" in
#' 82 bp of an 84% AT region could be a small indel misaligned as
#' substitutions. Reads are anchored on a unique, non-variable 20-mer just
#' right of the locus (CONTIG18:373,247-373,266; the left side is
#' repetitive AT sequence) and the 130 bp before the anchor in each read
#' (373,117-373,246 if no indel) are taken, on the reference strand.
#' Distinct sequences are counted per pool and compared with the assembly:
#' length difference (indel) and mismatches.
#'
#' Pools: three where mapped B is common (SN-2-97, SN-2-02, SN-1-97) and
#' three where it's rare (SN-1-99, SN-2-15, SN-2-08). Reads:
#' `/Volumes/lan_images4/ill/dna/trimmed/`.
#'
#' Writes `_data/contig18-flank.txt`.
#'

library(tidyverse)
library(seqinr)

out_file <- "_data/contig18-flank.txt"
sink(out_file, split = TRUE)

ctg <- "CONTIG18"
anc_start <- 373247
anc_len <- 20L
up_len <- 130L
trim_dir <- "/Volumes/lan_images4/ill/dna/trimmed"
pools <- c("SN-2-97_S37", "SN-2-02_S30", "SN-1-97_S61", "SN-1-99_S28",
           "SN-2-15_S34", "SN-2-08_S32")

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
cat("Assembly, 130 bp before the anchor:\n", ref_up, "\n", sep = "")
stopifnot(n_occ == 1)

#' Reads with the anchor (or its reverse complement), oriented to the
#' reference strand; the 130 bp before the anchor (reads that don't reach
#' that far are dropped).
pat <- tempfile(fileext = ".txt")
write_lines(c(anchor, revcomp(anchor)), pat)
reads <- map(pools, \(p) {
    tar <- file.path(trim_dir, paste0("trimmed_", p, ".tar"))
    members <- sprintf("trimmed_%s_L002_R%i_001.fastq.gz", p, 1:2)
    cmd <- sprintf("tar -xOf '%s' %s | gzip -dc | awk 'NR %% 4 == 2' | grep -F -f '%s'",
                   tar, paste(members, collapse = " "), pat)
    r <- system(cmd, intern = TRUE)
    r <- ifelse(grepl(anchor, r, fixed = TRUE), r, map_chr(r, revcomp))
    pos <- regexpr(anchor, r, fixed = TRUE)
    up <- ifelse(pos > up_len, substr(r, pos - up_len, pos - 1), NA_character_)
    tibble(pool = p, n_anchor = length(r), up = up[!is.na(up)])
}) |>
    list_rbind()

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
cat("\n== Most common sequences (all pools), 130 bp before the anchor\n")
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
    left_join(reads |> distinct(pool, n_anchor), by = "pool") |>
    arrange(match(pool, pools)) |>
    as.data.frame() |>
    print()

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
