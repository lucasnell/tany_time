# `01-variant-calling/`

Pool-seq reads → filtered, realigned BAMs → mpileup files.
One job per pool (`../samples/all-pools.txt`) except `05-kraken`. The
combined reference and the SILVA database come from `../00-references/`.
`../check-step.sh 01-variant-calling/<step>` summarizes any step (jobs, run
times, missing outputs, and for `02-bwa`, reads per reference group).

Reads are mapped competitively (as in DEST) to a combined reference: the
*T. gracilentus* assembly plus other species whose DNA could be in the
samples. Reads that match another species better map there, reads that
match both equally get low MAPQ and are filtered, and the mpileup includes
only *T. gracilentus* contigs.

| Step | Does | Input (`/staging/l/lnell/`) | Output (`/staging/l/lnell/dna/`) |
|---|---|---|---|
| `01-fastp` | Adapter and polyG trimming, overlap base correction (`--correction`), no quality filtering (bwa soft-clips low-quality ends); FastQC on trimmed reads | `dna/raw_fq/` | `trimmed/trimmed_<pool>.tar` (reads), `trimmed/trimmed_<pool>.tar.gz` (reports) |
| `02-bwa` | Merge overlapping pairs (`bbmerge.sh`); align merged and unmerged reads separately (`bwa mem -M`) to the combined reference; keep MAPQ ≥ 20 and primary alignments (`-q 20 -F 0x100`); merge (`gatk MergeSamFiles`); count reads per species | `trimmed/`, `hologenome.tar` | `bwa/<pool>_bwa.bam` (all species); logs, `bamtools stats`, `<pool>_species-counts.tsv`, and unmapped reads in `bwa/<pool>_bwa.tar.gz` |
| `03-mpileup` | Remove duplicates (`gatk MarkDuplicates`); realign around indels (GATK 3 `RealignerTargetCreator` + `IndelRealigner`, 9+ h), both on the combined reference; `samtools mpileup -B -Q 25` on *T. gracilentus* contigs only | `bwa/`, `hologenome.tar` | `mpileup/<pool>_mpileup.txt.gz`, `mpileup/<pool>_mpileup.tar.gz` |
| `04-bam-cov` | Per-base coverage (`bedtools genomecov -bga`) on *T. gracilentus* contigs, used in `_R/sample-qc.R` | `bwa/`, `hologenome.tar` | `bwa/<pool>_bwa_cov.bed.gz` |
| `05-kraken` | Diagnostic, not used for pool inclusion. Kraken2 (`kraken-env`; confidence 0.1) on each pool's unmapped reads: PlusPF-16 on a random 1M (overall composition, endosymbiotic bacteria, viruses) and SILVA on all of them (eukaryotic parasites: fungi, oomycetes, microsporidia, mermithids, mites, ...). Summarized by `_R/kraken-summary.R` | `bwa/`, `kraken2-db/` (PlusPF-16 from https://benlangmead.github.io/aws-indexes/k2; SILVA from `../00-references/03-silva-db`) | `kraken/kraken.tar` |

## Notes

- Reads per reference group (`<pool>_species-counts.tsv` from `02-bwa`)
  give a mitochondrial contamination index per trap taxon
  (`_R/mito-contam.R`; relative, because mitochondrial copy number varies
  among species).
- `samtools mpileup -Q 25` matches the base-quality threshold that
  `Mpileup2Sync.py` uses for counts in `02-allele-freqs/01-snape` (as in
  DESTv3). Outputs made before 2026-10-09 used the default (13).
- Before 2026-10-10, `03-mpileup` also wrote 500-bp window coverage
  (`*_mpileup_win500.txt.gz`); it was removed because overlapping windows
  undercounted and nothing used the output.
