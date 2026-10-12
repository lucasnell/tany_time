# `samples/`: job lists

Written by `_R/make-job-lists.R` (run from the repo root). Do not edit by
hand; edit the R script and re-run it.

| File | Contents | Used by |
|---|---|---|
| `all-pools.txt` | All 62 sequenced pools | `01-variant-calling/` (`01-fastp`, `02-bwa`, `03-mpileup`, `05-bam-cov`, `06-kraken`) |
| `all-pools-n-adults.txt` | `pool,n_adults` for all 62 pools (from `_data/full-sample-info.csv`) | `02-allele-freqs/` (`01-snape`, `02-npstat`) |
| `time-pools-candidates.txt` | All 30 time-series (SN) pools | `02-allele-freqs/03-pool-qc` |
| `space-pools-candidates.txt` | All 18 spatial (2019) pools | `02-allele-freqs/03-pool-qc` |
| `time-pools-cov30.txt`, `space-pools-cov30.txt` | Pools with median coverage ≥ 30 (the 2025 rule; 17 and 18) | sensitivity analyses |
| `time-pools.txt`, `space-pools.txt` | Pools included after QC; **written by `_R/pool-qc.R`**, not by `make-job-lists.R` (absent until QC rules are set) | `02-allele-freqs/04-combine-syncs` |

SNAPE and npstat run on all 62 pools, so changing which pools are included
only requires re-running `02-allele-freqs/04-combine-syncs` onward.

`n_adults` sets the pool size (`2 × n_adults` haploid genomes) for SNAPE,
npstat, and BayPass, so it must match `_data/full-sample-info.csv`.
