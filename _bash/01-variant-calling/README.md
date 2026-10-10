# `01-variant-calling/`

Pool-seq reads → filtered, realigned BAMs → mpileup files.
One job per pool (`../samples/all-pools.txt`) except `04-ont-pi`.

| Step | Does | Input (`/staging/lnell/`) | Output (`/staging/lnell/dna/`) |
|---|---|---|---|
| `01-fastp` | Adapter and polyG trimming, overlap base correction (`--correction`), no quality filtering (bwa soft-clips low-quality ends); FastQC on trimmed reads | `dna/raw_fq/` | `trimmed/trimmed_<pool>.tar` (reads), `trimmed/trimmed_<pool>.tar.gz` (reports) |
| `02-bwa` | Merge overlapping pairs (`bbmerge.sh`); align merged and unmerged reads separately (`bwa mem -M`); keep MAPQ ≥ 20 and primary alignments (`-q 20 -F 0x100`); merge (`gatk MergeSamFiles`) | `trimmed/`, `Tgraci_assembly.fasta.gz` | `bwa/<pool>_bwa.bam`; logs, `bamtools stats`, and unmapped reads in `bwa/<pool>_bwa.tar.gz` |
| `03-mpileup` | Remove duplicates (`gatk MarkDuplicates`); realign around indels (GATK 3 `RealignerTargetCreator` + `IndelRealigner`, 9+ h); `samtools mpileup -B -Q 25` | `bwa/`, assembly | `mpileup/<pool>_mpileup.txt.gz`, `mpileup/<pool>_mpileup.tar.gz` |
| `04-ont-pi` | Map Nanopore reads from one individual (minimap2), call variants (bcftools), mean π in 1-kb windows (vcftools). Result (0.0088) is the SNAPE prior θ | `dna/ont/basecalls_guppy-5.0.11.fastq.gz`, assembly | `ont/ont_variants.bcf`, `ont/ont_variants.windowed.pi.gz` |
| `05-bam-cov` | Per-base coverage (`bedtools genomecov -bga`), used in `_R/sample-qc.R` | `bwa/` | `bwa/<pool>_bwa_cov.bed.gz` |

## Notes

- `samtools mpileup -Q 25` matches the base-quality threshold that
  `Mpileup2Sync.py` uses for counts in `02-allele-freqs/01-snape` (as in
  DESTv3). Outputs made before 2026-10-09 used the default (13).
- Before 2026-10-10, `03-mpileup` also wrote 500-bp window coverage
  (`*_mpileup_win500.txt.gz`); it was removed because overlapping windows
  undercounted and nothing used the output.
