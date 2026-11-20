#'
#' What is in the genome at the CONTIG18 indel haplotype? (Owner,
#' 2026-11-20; exploratory.) The locus (CONTIG18:373,124-373,206; 9 mapped
#' "SNPs") is a haplotype, B, that differs from the assembly (A) by small
#' indels and substitutions (`_data/contig18-flank.txt`).
#'
#'   1. Gene context: AUGUSTUS features (`~/_data/midgenomes/features/
#'      Tgraci_features.gff3.gz`) and annotated repeats (RepeatMasker,
#'      `~/_data/midgenomes/repeats/Tgraci_repeats/Tgraci_repeats_locs.gff3`)
#'      within 1 kb; the locus's position relative to exons and splice
#'      sites.
#'   2. Sequence of A around the locus, composition, and homopolymer runs.
#'   3. B's sequence from the reads: reads anchored just right of the locus
#'      (cached by `_R/contig18-flank.R`) end at the same point, so they
#'      line up from the anchor backward; per position, the majority base
#'      over B reads (and over A reads, as a check against the assembly),
#'      where >= 5 reads reach.
#'   4. A vs. B alignment (blastn, both sequences ending at the anchor).
#'
#' Writes `_data/contig18-indel.txt`.
#'

library(tidyverse)
library(seqinr)

out_file <- "_data/contig18-indel.txt"
sink(out_file, split = TRUE)

ctg <- "CONTIG18"
lo <- 373124
hi <- 373206
anc_start <- 373247
anc_len <- 20L



# =============================================================================*
# ---- 1. gene and repeat context ----
# =============================================================================*


feat <- read_tsv("~/_data/midgenomes/features/Tgraci_features.gff3.gz",
                 comment = "#", col_types = "ccciicccc",
                 col_names = c("contig", "source", "type", "start", "end",
                               "score", "strand", "phase", "attr")) |>
    filter(contig == ctg, end >= lo - 1000, start <= hi + 1000)
cat("== 1. Features within 1 kb of", ctg, lo, "-", hi, "\n")
feat |>
    select(type, start, end, strand, attr) |>
    mutate(attr = str_trunc(attr, 60)) |>
    as.data.frame() |>
    print()
ex <- feat |> filter(type == "exon") |> arrange(start)
intr <- tibble(start = head(ex$end, -1) + 1, end = tail(ex$start, -1) - 1)
cat("\nIntron(s) within 1 kb (g4149 has 7 exons; these are exons 6 and 7):\n")
print(as.data.frame(intr |> mutate(length = end - start + 1)))
cat("Locus vs. the next exon's start (", ex$start[ex$start > hi][1],
    "): ", ex$start[ex$start > hi][1] - hi, " to ",
    ex$start[ex$start > hi][1] - lo, " bp upstream of it\n", sep = "")

rep <- read_tsv("~/_data/midgenomes/repeats/Tgraci_repeats/Tgraci_repeats_locs.gff3",
                comment = "#", col_types = "cccddcccc",
                col_names = c("contig", "source", "type", "start", "end",
                              "score", "strand", "phase", "attr")) |>
    filter(contig == ctg, end >= lo - 1000, start <= hi + 1000)
cat("\nAnnotated repeats within 1 kb:\n")
print(as.data.frame(rep |> select(start, end, strand, attr)))



# =============================================================================*
# ---- 2. assembly sequence ----
# =============================================================================*


asm <- read.fasta("~/_data/midgenomes/assemblies/Tgraci_assembly.fasta.gz",
                  seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
s <- toupper(asm[[ctg]][[1]])
rm(asm)
w_lo <- 373060
w_hi <- 373300
seq_a <- substr(s, w_lo, w_hi)
cat("\n== 2. Assembly (A), ", ctg, ":", w_lo, "-", w_hi, " (lower case: the",
    " locus's 9 mapped SNPs span ", lo, "-", hi, "):\n", sep = "")
v <- strsplit(seq_a, "")[[1]]
pos <- w_lo:w_hi
v[pos >= lo & pos <= hi] <- tolower(v[pos >= lo & pos <= hi])
cat(strwrap(paste(v, collapse = ""), 60), sep = "\n")
loc_seq <- substr(s, lo, hi)
cat("AT content: locus ", round(mean(strsplit(loc_seq, "")[[1]] %in%
                                       c("A", "T")), 3),
    "; 241-bp window ", round(mean(strsplit(seq_a, "")[[1]] %in%
                                     c("A", "T")), 3), "\n", sep = "")
runs <- rle(strsplit(substr(s, lo - 20, hi + 20), "")[[1]])
cat("Homopolymer runs >= 5 bp within 20 bp of the locus: ",
    paste0(runs$values, runs$lengths)[runs$lengths >= 5] |>
        paste(collapse = ", "), "\n", sep = "")



# =============================================================================*
# ---- 3. haplotype sequences from reads ----
# =============================================================================*


anchor <- substr(s, anc_start, anc_start + anc_len - 1)
raw <- read_rds("_data/contig18-flank-reads.rds")
up <- map_chr(raw$read, \(r) {
    p <- regexpr(anchor, r, fixed = TRUE)
    if (p <= 1) return(NA_character_)
    substr(r, 1, p - 1)
})
up <- up[!is.na(up) & nchar(up) >= 90]
ref90 <- substr(s, anc_start - 90, anc_start - 1)
b90 <- "ACTTTTGAATTTTTTAAATTTTTTTAAGAATTTATTTACGTTTTTAAATTTTTTTGAATTTTAAATTAATTTTTTTATTTTTTAGACACT"
last90 <- substr(up, nchar(up) - 89, nchar(up))
dA <- as.vector(adist(last90, ref90))
dB <- as.vector(adist(last90, b90))
hap <- ifelse(pmin(dA, dB) > 10, "other", ifelse(dA <= dB, "A", "B"))
cat("\n== 3. Reads anchored with >= 90 bp before the anchor:",
    length(up), "(A", sum(hap == "A"), ", B", sum(hap == "B"), ", other",
    sum(hap == "other"), ")\n")

#' Majority base by distance before the anchor, where >= 5 reads reach.
consensus <- function(x, min_n = 5) {
    rv <- map(x, \(r) rev(strsplit(r, "")[[1]]))
    L <- max(lengths(rv))
    out <- character(0)
    for (i in seq_len(L)) {
        b <- map_chr(rv, \(r) if (length(r) >= i) r[i] else NA_character_)
        b <- b[!is.na(b)]
        if (length(b) < min_n) break
        out <- c(out, names(sort(table(b), decreasing = TRUE))[1])
    }
    paste(rev(out), collapse = "")
}
cons_a <- consensus(up[hap == "A"])
cons_b <- consensus(up[hap == "B"])
ref_same <- substr(s, anc_start - nchar(cons_a), anc_start - 1)
cat("A consensus (", nchar(cons_a), " bp before the anchor) matches the ",
    "assembly: ", cons_a == ref_same, "\n", sep = "")
cat("B consensus (", nchar(cons_b), " bp before the anchor):\n", sep = "")
cat(strwrap(cons_b, 60), sep = "\n")



# =============================================================================*
# ---- 4. A vs. B alignment ----
# =============================================================================*


#' A over the same span, plus 20 bp on the left (B may be longer or
#' shorter than A over that span).
a_cmp <- substr(s, anc_start - nchar(cons_b) - 20, anc_start - 1)
fa <- tempfile(fileext = ".fa")
fb <- tempfile(fileext = ".fa")
write.fasta(list(a_cmp), "A_assembly", fa)
write.fasta(list(cons_b), "B_reads", fb)
al <- system2("blastn", c("-query", fb, "-subject", fa, "-dust", "no",
                          "-word_size", "7", "-outfmt", "0", "-evalue", "10"),
              stdout = TRUE)
cat("\n== 4. B (query) vs. A (subject, ", ctg, ":",
    anc_start - nchar(cons_b) - 20, "-", anc_start - 1, "); blastn\n",
    sep = "")
i <- grep("^ Score =", al)[1]
j <- grep("^ Score =", al)[2]
cat(al[i:(j - 1)], sep = "\n")
tab <- system2("blastn", c("-query", fb, "-subject", fa, "-dust", "no",
                           "-word_size", "7", "-evalue", "10",
                           "-outfmt", shQuote("6 qstart qend sstart send length mismatch gaps pident")),
               stdout = TRUE)
cat("\nBest hit (query start/end, subject start/end, length, mismatches,",
    "gaps, % identity):\n")
cat(tab[1], sep = "\n")
cat("B is", nchar(cons_b), "bp over the span that is",
    as.integer(str_split_1(tab[1], "\t")[4]) -
        as.integer(str_split_1(tab[1], "\t")[3]) + 1, "bp in A\n")

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
