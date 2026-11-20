#'
#' How far does the CONTIG18 divergent haplotype (B) extend? (Owner,
#' 2026-11-20; exploratory.) The nanopore reads come from one individual;
#' their variant calls (`00-references/04-ont-pi`: Dorado 0.9.6 sup reads
#' mapped to the assembly, bcftools SNP calls; copied from staging to
#' `tany_dir("ont")`) show that individual heterozygous at the locus's SNPs
#' with the B alleles (`_data/contig18-indel.txt`), i.e. it carries A and
#' B. So the density of its heterozygous sites along the contig traces
#' where A and B differ: a divergent stretch shows as a run of dense
#' heterozygous sites, against the genome-wide background.
#'
#' Heterozygous SNP calls (GT 0/1, QUAL >= 30) per 250-bp window,
#' CONTIG18:355,000-395,000 and genome-wide (for the background
#' distribution), and the calls within 2 kb of the locus with their
#' position relative to g4149's exons.
#'
#' Writes `_data/contig18-ont.txt`.
#'

library(tidyverse)

source("_R/data-dir.R")

out_file <- "_data/contig18-ont.txt"
sink(out_file, split = TRUE)

ctg <- "CONTIG18"
lo <- 373124
hi <- 373206
win <- 250

vcf <- read_tsv(tany_dir("ont", "ont_variants_dorado-0.9.6_sup.vcf.gz"),
                comment = "#", col_types = "cicccdcccc",
                col_names = c("contig", "pos", "id", "ref", "alt", "qual",
                              "filter", "info", "format", "sample")) |>
    mutate(gt = str_extract(sample, "^[^:]+"),
           indel = str_detect(info, "INDEL"))
het <- vcf |> filter(gt == "0/1", qual >= 30, !indel)
cat("Heterozygous SNP calls (QUAL >= 30): ", nrow(het), " genome-wide; ",
    "genome-wide density ", signif(nrow(het) / 91518485 * 1000, 3),
    " per kb\n", sep = "")

#' Background: heterozygous calls per 250-bp window over all windows of
#' all contigs (windows without calls count as 0).
lens <- read_lines("~/_data/midgenomes/assemblies/Tgraci_assembly.fasta.gz",
                   n_max = -1)
hdr <- which(startsWith(lens, ">"))
clen <- tibble(contig = str_remove(lens[hdr], "^>") |> word(1),
               len = c(diff(hdr) - 1, length(lens) - tail(hdr, 1)))
seq_len_by <- map_int(seq_along(hdr), \(i) {
    a <- hdr[i] + 1
    b <- if (i < length(hdr)) hdr[i + 1] - 1 else length(lens)
    sum(nchar(lens[a:b]))
})
clen$len <- seq_len_by
n_win <- sum(clen$len %/% win)
counts <- het |>
    mutate(w = pos %/% win) |>
    count(contig, w)
bg <- c(counts$n, rep(0L, n_win - nrow(counts)))
cat("Calls per ", win, "-bp window, genome-wide: mean ", signif(mean(bg), 3),
    "; 99th percentile ", quantile(bg, 0.99), "; 99.9th ",
    quantile(bg, 0.999), "; max ", max(bg), "\n", sep = "")

cat("\n== Calls per", win, "bp window,", ctg, "355-395 kb (start: count)\n")
prof <- tibble(w = (355000 %/% win):(395000 %/% win)) |>
    left_join(counts |> filter(contig == ctg), by = "w") |>
    mutate(n = replace_na(n, 0L), start = w * win,
           pct_rank = map_dbl(n, \(x) mean(bg < x)))
prof |>
    filter(n > 0) |>
    transmute(start, n, genome_percentile = round(100 * pct_rank, 1)) |>
    as.data.frame() |>
    print()

cat("\n== Heterozygous calls within 2 kb of the locus\n")
feat <- read_tsv("~/_data/midgenomes/features/Tgraci_features.gff3.gz",
                 comment = "#", col_types = "ccciicccc",
                 col_names = c("contig", "source", "type", "start", "end",
                               "score", "strand", "phase", "attr")) |>
    filter(contig == ctg, type == "exon", str_detect(attr, "g4149\\.t1")) |>
    arrange(start) |>
    mutate(exon = row_number())
cat("(g4149 exons:", paste(paste0(feat$start, "-", feat$end), collapse = ", "),
    ")\n")
het |>
    filter(contig == ctg, pos >= lo - 2000, pos <= hi + 2000) |>
    mutate(where = map_chr(pos, \(p) {
        e <- which(feat$start <= p & feat$end >= p)
        if (length(e) > 0) return(paste0("exon ", feat$exon[e]))
        i <- which(feat$end[-nrow(feat)] < p & feat$start[-1] > p)
        if (length(i) > 0) return(sprintf("intron %i-%i", i, i + 1))
        "outside g4149"
    }),
    in_locus_span = pos >= lo & pos <= hi,
    ad = str_extract(sample, "[^:]+$")) |>
    select(pos, ref, alt, qual, ad, where, in_locus_span) |>
    mutate(qual = round(qual)) |>
    as.data.frame() |>
    print()

#' Heterozygous-call density by segment of g4149 around the locus:
seg <- tibble(segment = c("exon 6", "intron 6-7 (outside the locus)",
                          "locus span", "exon 7"),
              start = c(feat$start[6], feat$end[6] + 1, lo, feat$start[7]),
              end = c(feat$end[6], feat$start[7] - 1, hi, feat$end[7])) |>
    mutate(n = map2_int(start, end, \(a, b) {
        sum(het$contig == ctg & het$pos >= a & het$pos <= b)
    }),
    length = end - start + 1,
    n = ifelse(segment == "intron 6-7 (outside the locus)",
               n - sum(het$contig == ctg & het$pos >= lo & het$pos <= hi), n),
    length = ifelse(segment == "intron 6-7 (outside the locus)",
                    length - (hi - lo + 1), length),
    per_kb = round(1000 * n / length, 1))
cat("\n== Heterozygous calls per kb by segment (genome-wide",
    signif(nrow(het) / 91518485 * 1000, 3), "per kb)\n")
print(as.data.frame(seg))

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
