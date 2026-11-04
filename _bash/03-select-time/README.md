# `03-select-time/`: BayPass

Tests for SNPs whose allele frequencies covary with *T. gracilentus*
abundance across the 18 time-series pools (`../samples/time-pools.txt`). Runs in `baypass-env`
(BayPass 3.1 in the v2.2.0 container; existing results used 2.31).

Everything runs for a grid of missing-data settings: the minimum proportion
of pools with data at a SNP (`MIN_PROP_POOLS` = 0.5, 0.75, 1.0; tags
`mp050`, `mp075`, `mp100`). All are reported; `mp075` is the main analysis
(decided 2026-10-13). Missingness tracks pool coverage, so it isn't random
with respect to abundance, and published thresholds vary widely.

| Step | Does | Jobs | Output (`/staging/l/lnell/dna/baypass/`) |
|---|---|---|---|
| `00-baypass-inputs` | Read `count_masked_time.sync.gz` with poolfstat (`popsync2pooldata`, default filters), with masked pools as 0 counts (BayPass treats these as missing data, manual v3.1 §2.2.1); keep SNPs with data in ≥ `MIN_PROP_POOLS` of pools and mean reference-allele read proportion (over pools with data) in (0.1, 0.9); write genotype, pool-size, SNP-detail, and logN covariate files, for all SNPs and as 3 thinned subsamples | one per setting (3; `baypass_inputs_info.txt`) | `baypass-inputs_<tag>.tar.gz`, `baypass-inputs_<tag>_3.tar.gz` |
| `01-baypass-init` | Core model to estimate Ω (population covariance matrix) on each of the 3 thinned subsamples (`sub1`–`sub3`) of each setting | 9 (`baypass_init_info.txt`) | `init/baypass-init-<tag>-sub<k>.tar.gz` |
| `02-check-omegas` | FMD distances (Förstner and Moonen 2003) between the subsamples' Ω within each setting; reports PASS if all pairs are < 0.5 (Gautier 2015), else WARN. Diagnostic only | 1 | `init/omega-check.tsv` |
| `03-baypass-beta` | Covariate model (`-auxmodel`) on each subsample with that subsample's Ω, for each Ising β (0, 1) and run (1–10); keeps the subsample's `snpdet` file for combining | 180 (`baypass_beta_info.txt`) | `beta/baypass-<tag>-sub<k>-beta<β>-run<r>.tar` |

Subsamples (owner, 2026-11-04): `00-baypass-inputs` deals SNPs in genome
order to 3 subsamples in turn (poolfstat thinning), so they don't overlap
and together hold every SNP. A full-set core-model run took ~71 h for
`mp050` (1.02 million SNPs), at CHTC's 72-h limit. Results are combined
across subsamples locally. Neighboring SNPs in a subsample are ~270 bp
apart (`mp050`), so subsamples still carry LD.

## Inputs

- `00-baypass-inputs` needs `/staging/l/lnell/tany-sample-info.tar.gz`,
  containing `_data/full-sample-info.csv` and `_data/tany-abundances.csv`
  (see comments at the top of the script for how to make and upload it), and
  `snape/count_masked_time.tar` from `02-allele-freqs/04-combine-syncs`.
- Job lists `baypass_inputs_info.txt` (`tag,min_prop_pools`),
  `baypass_init_info.txt` (`tag,subsample`), and `baypass_beta_info.txt`
  (`ising_beta,run,seed,tag,subsample`) are written by
  `_R/make-job-lists.R`. Seeds reproduce the 2025 runs and are shared
  across settings and subsamples.
- Seeds for `01-baypass-init` are fixed in the script.

## Before re-running (see `TODO.md`)

- Both BayPass scripts set `-nval` (`NVAL=5000`; 3.1's default is 250).
  Expect roughly 25–50 h per job with 8 CPUs (2.31 runs with ≤ 35,000
  iterations took 9–12 h; 5,000 samples under 3.1 defaults is ~102,000
  iterations). That's 6 Ω jobs plus 60 covariate jobs. CHTC's default run
  limit is 72 h. Measured 2026-11-04: ~2.5 s per iteration on 1.02 million
  SNPs (~71 h), hence the subsamples. Check the first covariate job's
  progress before submitting the rest.
  BayPass 3.1 scales covariates by default (2.31 runs used raw logN).
- Existing Ω came from 2.31 with 8 threads, affected by a bug fixed in 2.4.
- Local analysis of outputs: `_R/baypass-check-omegas.R`,
  `_R/baypass-analysis.R`, `_R/baypass-check-candidates.R`. These read the
  2025 output names (no setting tag) and need updating for the grid.
