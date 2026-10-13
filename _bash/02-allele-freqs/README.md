# `02-allele-freqs/`

mpileup files → SNAPE-pooled allele frequencies, masked allele counts, and
diversity statistics. Follows the DEST pipeline (Kapun et al. 2021,
doi:10.1093/molbev/msab259) using its helper scripts (in the container's
`/app/bin`).

| Step | Does | Jobs | Output (`/staging/l/lnell/dna/`) |
|---|---|---|---|
| `01-snape` | Uppercase reference and pileup; counts, coverage, and indel positions (`Mpileup2Sync.py`, base quality ≥ 25); SNAPE-pooled per contig; convert to gSYNC (`SNAPE2SYNC.py`); mask (`MaskSYNC_snape.py`) twice (below) | one per pool (`../samples/all-pools-n-adults.txt`) | `snape/<pool>_snape.tar` |
| `02-npstat` | Convert the partially masked gSYNC to pileup; npstat per contig in 10-kb windows (π, Watterson's θ, Tajima's D) at sites in the SNP file, with `-maxcov` set to each contig's maximum coverage (npstat's setup grows with `maxcov`²; output unchanged, `../tests/npstat-maxcov.sh`) | one per pool | `npstat/<pool>_npstat_10kb.stat.gz` |
| `03-pool-qc` | Per-pool QC metrics for all candidate pools (below) | one | `pool-qc/pool-qc.tar` |
| `04-combine-syncs` | Combine partially masked gSYNC files across included pools for time, space, and all pools, keeping sites SNAPE calls polymorphic (P ≥ 0.9) in at least one pool | one | `snape/{snape,count}_masked_{time,space,all}.tar` |

## SNAPE settings (`01-snape`)

- `-nchr` = 2 × `n_adults`; θ = 0.0088 (from `01-variant-calling/04-ont-pi`);
  D = 0.01; informative prior; unfolded.
- Masking (`MaskSYNC_snape.py`) removes sites with coverage < 10 or above
  the 95th percentile of the contig, within repeats
  (`Tgraci_repeats_locs.gff3.gz`), or within 5 bp of an indel
  (as in Kapun et al. 2021).
- **Fully masked** (`<pool>_snape_masked.sync`): also removes sites with
  P(polymorphic) < 0.9. This differs from DESTv3, which removes only
  0.1 < P < 0.9 (i.e., it keeps confidently monomorphic sites); the local
  version follows the GonzalezLab SNP-caller benchmark (see comments in
  `MaskSYNC_snape.py`). Positions kept here form the npstat SNP file.
- **Partially masked** (`<pool>_snape_part_masked.sync`): no P(polymorphic)
  filter; input to npstat (so monomorphic sites count toward diversity) and
  to `04-combine-syncs`.

## Choosing pools (`03-pool-qc`, then `_R/pool-qc.R`)

Pools are included or excluded from QC metrics, decided before any
allele-frequency analysis (following the DEST approach: Kapun et al. 2021;
Nunez et al. 2025). `03-pool-qc` computes, for the candidate pools in
`../samples/{time,space}-pools-candidates.txt`:

| Metric | Source |
|---|---|
| Missing-data rate after masking | partially masked gSYNC (`pool-qc-sync.py`) |
| Number of SNPs and private SNPs (P ≥ 0.9 only in that pool, among the set's candidates) | partially masked gSYNC (`pool-qc-sync.py`) |
| Duplication rate | Picard MarkDuplicates report (`03-mpileup`); fastp estimate also reported |
| GC content | fastp (`01-fastp`) |
| Mapping rate to *T. gracilentus* (`target_rate`: reads assigned to its contigs, MAPQ ≥ 20, primary / all records); rate to other species, any-genome mapping rate, and usable rate are also reported | `bamtools stats` and species counts (`02-bwa`) |

pN/pS isn't used because the *T. gracilentus* annotation isn't good enough.
`_R/pool-qc.R` excludes pools that are outliers in a robust PCA of the five
metrics (ROBPCA, with a threshold calibrated by parametric bootstrap for a
2.5% per-pool false-positive rate; see the script) or that lack a metric,
and writes `../samples/{time,space}-pools.txt`. The 2025 inclusion rule (median
coverage ≥ 30) is kept as `../samples/*-pools-cov30.txt` for sensitivity
analyses.

## Combined outputs (`04-combine-syncs`)

Inputs are the partially masked files. A site is kept if SNAPE calls it
polymorphic (P ≥ 0.9, `MIN_PP`) in at least one pool of the set, so pools
where an allele is fixed or lost still contribute their counts. Each tar
has `<name>.sync.gz` (sites passing that filter), `<name>.names.gz` (pool
order), and `<name>_noblanks.sync.gz` (those sites with data in every pool
after coverage, repeat, and indel masks).

- `snape_masked_*`: SNAPE info per pool (`ref count:alt count:ref quality:alt
  quality:P(polymorphic):P(fixed):mean minor-allele frequency`). Used by
  `_R/tidy-freq-data.R` and `_R/make-time-series-data.R`.
- `count_masked_*`: allele counts per pool (`A:T:C:G:N:del`). Where SNAPE
  calls a pool confidently monomorphic at a site (P(polymorphic) ≤ 0.1,
  `MAX_PP_MONO`), all of its reads are assigned to the reference base,
  removing sequencing-error reads (the reference-allele case of DESTv3's
  monomorphic filter). This matters because sites can now include pools
  fixed for the reference. True rare alleles that SNAPE calls monomorphic
  are erased too. The `time`
  file (with masked pools) is the BayPass input
  (`03-select-time/00-baypass-inputs`), where masked pools become missing
  data.

Before 2026-10-10 the inputs were the fully masked files, so `_noblanks`
kept only sites SNAPE called polymorphic in **every** pool, excluding sites
that drifted to fixation or near-loss in any one pool. All results made
before then (including BayPass candidates) used that SNP set.

Counts in these files are SNAPE's reference and most common alternative
allele counts per pool, so a pool's alternative allele can differ from
another pool's; poolfstat drops sites with more than two alleles overall.
