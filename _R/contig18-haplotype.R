#'
#' Where do reads carrying the CONTIG18 haplotype go? (Exploratory, owner,
#' 2026-11-19; see `_R/contig18-locus.R`.) Read depth at the locus falls as
#' the non-reference haplotype (B) gets more common, so B reads may be lost
#' at mapping. `02-bwa` maps to a hologenome (T. gracilentus plus
#' C. riparius, C. tentans, S. aterrima, P. akamusi) and keeps reads with
#' MAPQ >= 20, primary only: B reads would be lost if B matches another
#' place in the T. gracilentus assembly (a paralog) or another species'
#' genome about as well as CONTIG18.
#'
#' Builds the reference haplotype (A) and B (the 9 SNPs' B alleles;
#' `_data/contig18-locus.txt`) over CONTIG18:372,974-373,356 (150 bp
#' either side), and BLASTs both (blastn) against the T. gracilentus
#' assembly and every other assembly in `~/_data/midgenomes/assemblies/`.
#'
#' Writes `_data/contig18-haplotype.txt`.
#'

library(tidyverse)
library(seqinr)

source("_R/data-dir.R")

out_file <- "_data/contig18-haplotype.txt"
sink(out_file, split = TRUE)

ctg <- "CONTIG18"
lo <- 373124
hi <- 373206
w_lo <- lo - 150
w_hi <- hi + 150
asm_dir <- "~/_data/midgenomes/assemblies"
tmp <- tempfile("c18")
dir.create(tmp)

#' SNP alleles (mp050 subsamples) and B, as oriented in
#' `_R/contig18-locus.R`: B is the reference allele at 373,124 and the
#' alternate allele at the other 8 SNPs (poolfstat's labels).
#' BLAST runs with low-complexity filtering off (the locus is AT-rich).
inp <- tany_dir("baypass-inputs", "baypass-inputs_mp050_3.tar.gz")
untar(inp, files = paste0("baypass-inputs_mp050/tany.snpdet.sub", 1:3),
      exdir = tmp)
snps <- map(1:3, \(k) read_table(file.path(tmp, "baypass-inputs_mp050",
                                           paste0("tany.snpdet.sub", k)),
                                 col_names = c("contig", "pos", "ref", "alt"),
                                 col_types = "cdcc")) |>
    list_rbind() |>
    filter(contig == ctg, pos >= lo, pos <= hi) |>
    arrange(pos) |>
    mutate(b = ifelse(pos == lo, ref, alt))

asm <- read.fasta(file.path(asm_dir, "Tgraci_assembly.fasta.gz"),
                  seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
seq_a <- substr(asm[[ctg]][[1]], w_lo, w_hi)
seq_b <- strsplit(seq_a, "")[[1]]
snps$assembly <- NA_character_
for (i in seq_len(nrow(snps))) {
    j <- snps$pos[i] - w_lo + 1
    #' `ref`/`alt` here are poolfstat's labels, not necessarily the
    #' assembly's base, which must be one of them:
    stopifnot(toupper(seq_b[j]) %in% toupper(c(snps$ref[i], snps$alt[i])))
    snps$assembly[i] <- toupper(seq_b[j])
    seq_b[j] <- snps$b[i]
}
seq_b <- paste(seq_b, collapse = "")
rm(asm)
cat("\nAssembly base at each SNP (A = the assembly's haplotype):\n")
print(as.data.frame(snps |> mutate(b_is_assembly = b == assembly)))
cat("\nAssembly sequence, ", ctg, ":", w_lo, "-", w_hi, ":\n", sep = "")
cat(strwrap(seq_a, 80), sep = "\n")
cat("AT content: ", round(mean(strsplit(toupper(seq_a), "")[[1]] %in%
                               c("A", "T")), 3), "\n", sep = "")
q_file <- file.path(tmp, "query.fa")
write.fasta(list(seq_a, seq_b), c("A_reference", "B_haplotype"), q_file,
            nbchar = 80)
cat("\nQuery: ", ctg, ":", w_lo, "-", w_hi, " (", nchar(seq_a),
    " bp; A and B differ at ", nrow(snps), " sites)\n", sep = "")

#' BLAST against each assembly (database built in a temporary folder).
asms <- list.files(asm_dir, "_assembly\\.fasta\\.gz$", full.names = TRUE)
hits <- map(asms, \(f) {
    code <- str_remove(basename(f), "_assembly\\.fasta\\.gz$")
    fa <- file.path(tmp, paste0(code, ".fa"))
    system2("sh", c("-c", shQuote(paste("gzip -dc", shQuote(path.expand(f)),
                                        ">", shQuote(fa)))))
    system2("makeblastdb", c("-in", fa, "-dbtype", "nucl", "-out",
                             file.path(tmp, code)), stdout = FALSE)
    out <- system2("blastn", c("-query", q_file, "-db", file.path(tmp, code),
                               "-outfmt", shQuote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
                               "-evalue", "1e-10", "-max_target_seqs", "20",
                               "-dust", "no", "-soft_masking", "false"),
                   stdout = TRUE)
    file.remove(fa)
    if (length(out) == 0) return(NULL)
    read_tsv(I(out), col_names = c("query", "subject", "pident", "length",
                                   "mismatch", "gapopen", "qstart", "qend",
                                   "sstart", "send", "evalue", "bitscore"),
             col_types = "ccdiiiiiiidd") |>
        mutate(assembly = code, .before = 1)
}) |>
    list_rbind()

cat("\n== BLAST hits (e-value <= 1e-10), best 5 per query and assembly\n")
hits |>
    group_by(assembly, query) |>
    slice_max(bitscore, n = 5, with_ties = FALSE) |>
    ungroup() |>
    arrange(assembly, query, desc(bitscore)) |>
    as.data.frame() |>
    print()
cat("\nAssemblies searched:", paste(str_remove(basename(asms),
                                              "_assembly\\.fasta\\.gz$"),
                                   collapse = ", "), "\n")

unlink(tmp, recursive = TRUE)
sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
