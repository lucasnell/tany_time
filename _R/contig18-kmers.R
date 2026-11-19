#'
#' Superseded, not run to completion (2026-11-19): the flanking-read
#' sequences (`_R/contig18-flank.R`) show the B haplotype differs from the
#' assembly by small indels as well as substitutions, so the 25-mers built
#' here by substituting the 9 SNP alleles don't match B reads.
#'
#' Haplotype frequencies at the CONTIG18 locus straight from the reads,
#' without mapping (exploratory, owner, 2026-11-19; see
#' `_R/contig18-locus.R` and `_R/contig18-haplotype.R`). Mapped depth there
#' falls as the non-reference haplotype (B) gets more common, and B has no
#' match elsewhere in any assembly, so B reads may be lost to the MAPQ >= 20
#' filter in `02-bwa` (9 mismatches in 82 bp in an 84% AT region). Exact
#' 25-mers that span several SNPs are counted in each pool's trimmed reads
#' (both strands):
#'   k1  CONTIG18:373,155-373,179 (SNPs at 162, 167, 170)
#'   k2  CONTIG18:373,182-373,206 (SNPs at 184, 192, 196, 206)
#'   ctl CONTIG18:373,290-373,314 (no SNP; read coverage near the locus)
#' each as the assembly haplotype (A) and as B. Raw B frequency =
#' B / (A + B) per k-mer, compared with the mapped frequency from the
#' BayPass inputs (`_data/contig18-locus.txt`). Reads carrying a mix of A
#' and B alleles within a k-mer, or a sequencing error in it, aren't
#' counted.
#'
#' Pools: three where mapped B is common (SN-2-97, SN-2-02, SN-1-97), three
#' where it's rare (SN-1-99, SN-2-15, SN-2-08), and the two pools masked at
#' the locus (SN-1-95, SN-1-07). Reads: `/Volumes/lan_images4/ill/dna/trimmed/`.
#'
#' Writes `_data/contig18-kmers.txt`.
#'

library(tidyverse)
library(seqinr)

out_file <- "_data/contig18-kmers.txt"
sink(out_file, split = TRUE)

ctg <- "CONTIG18"
trim_dir <- "/Volumes/lan_images4/ill/dna/trimmed"
pools <- c("SN-2-97_S37", "SN-2-02_S30", "SN-1-97_S61", "SN-1-99_S28",
           "SN-2-15_S34", "SN-2-08_S32", "SN-1-95_S59", "SN-1-07_S23")
#' B alleles (`_data/contig18-haplotype.txt`): position -> base.
b_alleles <- c(`373124` = "A", `373132` = "T", `373162` = "C", `373167` = "G",
               `373170` = "T", `373184` = "G", `373192` = "T", `373196` = "G",
               `373206` = "T")
kmers_at <- tibble(name = c("k1", "k2", "ctl"),
                   start = c(373155, 373182, 373290)) |>
    mutate(end = start + 24)

asm <- read.fasta("~/_data/midgenomes/assemblies/Tgraci_assembly.fasta.gz",
                  seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
s <- toupper(asm[[ctg]][[1]])
rm(asm)
revcomp <- \(x) chartr("ACGT", "TGCA", paste(rev(strsplit(x, "")[[1]]),
                                             collapse = ""))
kmers <- kmers_at |>
    mutate(A = map2_chr(start, end, \(a, b) substr(s, a, b)),
           B = pmap_chr(list(start, end, A), \(a, b, x) {
               v <- strsplit(x, "")[[1]]
               p <- as.integer(names(b_alleles))
               i <- which(p >= a & p <= b)
               v[p[i] - a + 1] <- b_alleles[i]
               paste(v, collapse = "")
           })) |>
    pivot_longer(c(A, B), names_to = "hap", values_to = "seq") |>
    filter(!(name == "ctl" & hap == "B"))
cat("== 25-mers\n")
print(as.data.frame(kmers))

pat_file <- tempfile(fileext = ".txt")
pats <- kmers |>
    mutate(rc = map_chr(seq, revcomp)) |>
    pivot_longer(c(seq, rc), names_to = "strand", values_to = "kmer")
write_lines(unique(pats$kmer), pat_file)

counts <- map(pools, \(p) {
    tar <- file.path(trim_dir, paste0("trimmed_", p, ".tar"))
    members <- sprintf("trimmed_%s_L002_R%i_001.fastq.gz", p, 1:2)
    cmd <- sprintf(paste("tar -xOf '%s' %s | gzip -dc | awk 'NR %% 4 == 2' |",
                         "grep -F -o -f '%s' | sort | uniq -c"),
                   tar, paste(members, collapse = " "), pat_file)
    out <- system(cmd, intern = TRUE)
    x <- str_match(str_trim(out), "^(\\d+) ([ACGT]+)$")
    tibble(pool = p, kmer = x[, 3], n = as.integer(x[, 2]))
}) |>
    list_rbind()

tab <- pats |>
    left_join(counts, by = "kmer", relationship = "many-to-many") |>
    filter(!is.na(pool)) |>
    group_by(pool, name, hap) |>
    summarize(n = sum(n), .groups = "drop") |>
    complete(pool = pools, name, hap, fill = list(n = 0L)) |>
    filter(!(name == "ctl" & hap == "B"))

cat("\n== Read counts per 25-mer (both strands)\n")
tab |>
    pivot_wider(names_from = c(name, hap), values_from = n) |>
    mutate(raw_B_k1 = round(k1_B / (k1_A + k1_B), 3),
           raw_B_k2 = round(k2_B / (k2_A + k2_B), 3),
           k_total_over_ctl = round((k1_A + k1_B + k2_A + k2_B) / 2 / ctl_A, 2)) |>
    arrange(match(pool, pools)) |>
    as.data.frame() |>
    print()
cat("\n(k_total_over_ctl: A + B reads per locus 25-mer over reads at the",
    "control 25-mer; ~1 if no reads carry other variants there)\n")

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
