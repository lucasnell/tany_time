# `01-variant-calling/`

Pool-seq reads → filtered, realigned BAMs → mpileup files.
One job per pool (`../samples/all-pools.txt`) except `00-hologenome`,
`04-ont-pi`, and `06-kraken`.

Reads are mapped competitively (as in DEST) to a combined reference: the
*T. gracilentus* assembly plus other species whose DNA could be in the
samples. Reads that match another species better map there, reads that
match both equally get low MAPQ and are filtered, and the mpileup includes
only *T. gracilentus* contigs.

| Step | Does | Input (`/staging/lnell/`) | Output (`/staging/lnell/dna/`) |
|---|---|---|---|
| `00-hologenome` | Combine `Tgraci_assembly.fasta.gz` with other assemblies (currently *Chironomus riparius* and *C. tentans*, as proxies for the *Chironomus* species at Mývatn), renaming other sequences `nt_<code>__<name>`; bwa index, `.fai`, `.dict`, and a BED of target contigs | `Tgraci_assembly.fasta.gz`, `other-genomes/` | `/staging/lnell/hologenome.tar` |
| `01-fastp` | Adapter and polyG trimming, overlap base correction (`--correction`), no quality filtering (bwa soft-clips low-quality ends); FastQC on trimmed reads | `dna/raw_fq/` | `trimmed/trimmed_<pool>.tar` (reads), `trimmed/trimmed_<pool>.tar.gz` (reports) |
| `02-bwa` | Merge overlapping pairs (`bbmerge.sh`); align merged and unmerged reads separately (`bwa mem -M`) to the combined reference; keep MAPQ ≥ 20 and primary alignments (`-q 20 -F 0x100`); merge (`gatk MergeSamFiles`); count reads per species | `trimmed/`, `hologenome.tar` | `bwa/<pool>_bwa.bam` (all species); logs, `bamtools stats`, `<pool>_species-counts.tsv`, and unmapped reads in `bwa/<pool>_bwa.tar.gz` |
| `03-mpileup` | Remove duplicates (`gatk MarkDuplicates`); realign around indels (GATK 3 `RealignerTargetCreator` + `IndelRealigner`, 9+ h), both on the combined reference; `samtools mpileup -B -Q 25` on *T. gracilentus* contigs only | `bwa/`, `hologenome.tar` | `mpileup/<pool>_mpileup.txt.gz`, `mpileup/<pool>_mpileup.tar.gz` |
| `04-ont-pi` | Map Nanopore reads from one individual (minimap2), call variants (bcftools), mean π in 1-kb windows (vcftools). Result (0.0088) is the SNAPE prior θ | `dna/ont/basecalls_guppy-5.0.11.fastq.gz`, assembly | `ont/ont_variants.bcf`, `ont/ont_variants.windowed.pi.gz` |
| `05-bam-cov` | Per-base coverage (`bedtools genomecov -bga`) on *T. gracilentus* contigs, used in `_R/sample-qc.R` | `bwa/`, `hologenome.tar` | `bwa/<pool>_bwa_cov.bed.gz` |
| `06-kraken` | Diagnostic: Kraken2 on a random 1M unmapped reads per pool (`kraken-env`; confidence 0.1) to identify non-target DNA. Not used for pool inclusion | `bwa/`, Kraken2 database (`kraken2-db/`; PlusPF-16 from https://benlangmead.github.io/aws-indexes/k2) | `kraken/kraken.tar` |

## Notes

- Before running `00-hologenome`, copy the other assemblies
  (`~/_data/midgenomes/assemblies/{Cripar,Ctenta}_assembly.fasta.gz`) to
  `/staging/lnell/other-genomes/`. Add species to `OTHERS` in the script.

- `samtools mpileup -Q 25` matches the base-quality threshold that
  `Mpileup2Sync.py` uses for counts in `02-allele-freqs/01-snape` (as in
  DESTv3). Outputs made before 2026-10-09 used the default (13).
- Before 2026-10-10, `03-mpileup` also wrote 500-bp window coverage
  (`*_mpileup_win500.txt.gz`); it was removed because overlapping windows
  undercounted and nothing used the output.
