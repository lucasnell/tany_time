# `samples/`: job lists

Written by `_R/make-job-lists.R` (run from the repo root). Do not edit by
hand; edit the R script and re-run it.

| File | Contents | Used by |
|---|---|---|
| `all-pools.txt` | All 62 sequenced pools | `01-variant-calling/` (`01-fastp`, `02-bwa`, `03-mpileup`, `05-bam-cov`) |
| `all-pools-n-adults.txt` | `pool,n_adults` for all 62 pools (from `_data/full-sample-info.csv`) | `02-allele-freqs/` (`01-snape`, `02-npstat`) |
| `time-pools.txt` | 17 time-series pools passing QC | `02-allele-freqs/03-combine-syncs` |
| `space-pools.txt` | 18 spatial pools passing QC | `02-allele-freqs/03-combine-syncs` |

QC passing means median coverage ≥ 30 (`_R/sample-qc.R`). The reason for this
cutoff isn't recorded and it may change (see `TODO.md`); update the lists in
`_R/make-job-lists.R` if it does. SNAPE and npstat run on all 62 pools so that
changing the cutoff only requires re-running `03-combine-syncs` onward.

`n_adults` sets the pool size (`2 × n_adults` haploid genomes) for SNAPE,
npstat, and BayPass, so it must match `_data/full-sample-info.csv`.
