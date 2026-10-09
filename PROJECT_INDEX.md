# Project index

## Pipeline order

1. `_bash/01-variant-calling/` (cluster): `01-fastp` trim → `02-bwa` merge
   overlapping pairs with bbmerge, align with `bwa mem` → `03-mpileup` mark
   duplicates, GATK indel realignment, `samtools mpileup` → `04-ont-pi`
   (priors for SNAPE) → `05-bam-cov` (coverage `.bed.gz`).
2. `_R/tidy-sample-data.R` → `_data/full-sample-info.csv` (pool ID, date,
   lake, site, number of adults, lat/lon).
3. `_R/sample-qc.R`: picks pools with median coverage ≥ 30 (17 SN time pools,
   18 spatial pools; lists are in the script).
4. `_bash/02-allele-freqs/` (cluster): `01-snape` (SNAPE-pooled → masked
   gSYNC; undoes softmasking) → `02-npstat` → `03-combine-snapes` and
   `04-combine-counts` (combine passing pools into `time`, `space`, and `all`
   sets).
5. `_R/tidy-pop-abunds.R` → `_data/tany-abundances.csv` (year, generation,
   N, logN = log1p(N)).
6. `_bash/03-select-time/` (cluster):
   - `00-baypass-inputs`: poolfstat filtering, drops SNPs with mean
     reference-allele read proportion outside (0.1, 0.9), writes genotype,
     pool-size, and logN covariate files. 433,360 SNPs.
   - `01-baypass-init`: Ω (population covariance matrix) from all SNPs
     (`sub0`) and 3 subsamples.
   - `02-baypass-beta`: covariate model with Ω from all SNPs; 10 runs with
     different seeds; Ising β given as an argument. Current runs use
     `-auxmodel` with β = 0 (no dependency between neighbouring SNPs),
     because the assembly is not chromosome-level (contig-to-contig
     linkage is unknown). 1,000 retained MCMC samples per run.
7. `_R/baypass-check-omegas.R`: subsample Ω matrices disagree (FMD
   distance > 3.7 for subsample 3), so the all-SNP Ω is used.
8. `_R/baypass-analysis.R`: median Bayes factor (dB) across 10 `beta0` runs;
   candidates at ≥ 30 dB → `_data/candi_df.rds`; nearest gene within 5 kb →
   `_data/candi_funcs_df.rds`; BLASTP of candidate genes against
   *D. melanogaster* UniProt; GO over-representation with `clusterProfiler`.
9. `_R/make-time-series-data.R` → `_data/timeseries-data.csv` (candidate SNP
   frequencies × pools, with abundance).
10. `_R/make-time-series-plot.R`: frequency-vs-abundance plots,
    focal genes sox21b and SNF4Aγ → `_plots/focal-freqs.pdf`.

## `_R/`

| File | Role |
|---|---|
| `tidy-sample-data.R` | Builds `full-sample-info.csv` from archive and other-lakes sheets |
| `sample-qc.R` | MultiQC + coverage summaries; chooses pools |
| `tidy-pop-abunds.R` | Builds `tany-abundances.csv` |
| `tidy-freq-data.R` | Biallelic, ≥ 500 bp-thinned SNAPE frequencies with effective sample size → `allele-freqs-time.rds` |
| `snape-helpers.cpp` | Rcpp parsers for SNAPE/count sync strings, distance thinning |
| `crash-test.R` | Legacy: earlier binomial GLMM approach (lme4), superseded by BayPass |
| `baypass-format.R` | Local precursor to `00-baypass-inputs.sh` |
| `baypass-utils.R` | Copied BayPass utility functions |
| `baypass-check-omegas.R` | Ω consistency across subsamples |
| `baypass-analysis.R` | Candidate SNPs, gene matching, BLAST, enrichment |
| `baypass-check-candidates.R` | Confirms `candi_df.rds` matches local BayPass output; per-candidate run consistency → `_data/baypass-check-candidates.txt` |
| `make-time-series-data.R` | Builds `timeseries-data.csv` |
| `make-time-series-plot.R` | Candidate frequency vs. abundance plots |

## `_data/`

Tracked CSVs (column descriptions for the two sample sheets are in
`_data/README.md`):

- `tany-from-archive.csv`, `tany-other-lakes.csv`,
  `tany-other-lakes-gps.csv`, `trap_locations.csv`, `biotech_ids.csv` — raw
  sample records.
- `full-sample-info.csv` — combined sample table.
- `chironomid-abundances.csv` → `tany-abundances.csv` — emergence abundance
  by year and generation.
- `timeseries-data.csv` — 33 SNPs × 17 pools = 561 rows; columns described at
  the end of `make-time-series-data.R`.

Local-only (`.rds`, git-ignored): `bed-summs.rds`, `allele-freqs-time.rds`,
`allele-freqs-time-crash.rds`, `crash-mod{1,3}-cV.rds`, `candi_df.rds`,
`candi_funcs_df.rds`, `goall_map.rds`.

## Current numbers

From `_R/baypass-check-candidates.R` → `_data/baypass-check-candidates.txt`
(run 2026-10-09):

- 433,360 SNPs tested; 33 candidates (median ≥ 30 dB); saved
  `candi_df.rds` matches.
- Each candidate passes 30 dB in only 5–7 of 10 runs (12, 15, and 6 SNPs
  respectively).
- Maximum Bayes factor in any run is 52.96 dB, shared by the top
  candidates. This is a ceiling: it occurs when a SNP's posterior inclusion
  probability is 1.0 in all 1,000 retained MCMC samples
  (`tany-beta0-run1_summary_betai.out.gz`).

From `_data/candi_funcs_df.rds`: 26 of 33 candidates fall inside or within
5 kb of a gene (24 unique genes); 8 unique *D. melanogaster* UniProt best hits
at e < 1e-4 (listed with annotations in `baypass-analysis.R`).

## `_legacy/`

Older analyses moved from the repo root (2026-10-09); not maintained and not
runnable as-is (see `CLAUDE.md`). Includes sample summaries, spatial FST and
connectivity (`fst-dist.stan`), admixture, temporal diversity, diversity
simulations (`sim_midges.cpp` with the `pcg/` RNG headers), and genome-report
maps. Git-ignored outputs `fst-dist.rds` and `sample_summary.pdf` moved with
them.
