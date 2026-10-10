# `03-select-time/`: BayPass

Tests for SNPs whose allele frequencies covary with *T. gracilentus*
abundance across the 17 time-series pools. Runs in `baypass-env`
(BayPass 3.1 in the v2.0.0 container; existing results used 2.31).

| Step | Does | Jobs | Output (`/staging/lnell/dna/baypass/`) |
|---|---|---|---|
| `00-baypass-inputs` | Read `count_masked_time.sync.gz` with poolfstat (`popsync2pooldata`, default filters), with masked pools as 0 counts (BayPass treats these as missing data, manual v3.1 §2.2.1); keep SNPs with data in ≥ `MIN_PROP_POOLS` of pools (provisionally 0.75) and mean reference-allele read proportion (over pools with data) in (0.1, 0.9); write genotype, pool-size, SNP-detail, and logN covariate files, for all SNPs and as 3 thinned subsamples | one | `baypass-inputs.tar.gz`, `baypass-inputs_3.tar.gz` |
| `01-baypass-init` | Core model to estimate Ω (population covariance matrix): all SNPs (`sub0`) and each subsample (`sub1`–`sub3`) | four (`queue 4`) | `init/baypass-init-sub<k>.tar.gz` |
| `02-baypass-beta` | Covariate model (`-auxmodel`) with Ω from `sub0`, for each Ising β and run | one per line of `baypass_beta_info.txt` | `beta/baypass-beta<β>-run<r>.tar` |

## Inputs

- `00-baypass-inputs` needs `/staging/lnell/tany-sample-info.tar.gz`,
  containing `_data/full-sample-info.csv` and `_data/tany-abundances.csv`
  (see comments at the top of the script for how to make and upload it), and
  `snape/count_masked_time.tar` from `02-allele-freqs/04-combine-syncs`.
- `baypass_beta_info.txt` (`ising_beta,run,seed`; β ∈ {0, 1}, 10 runs each)
  is written by `_R/make-job-lists.R` and reproduces the 2025 seeds.
- Seeds for `01-baypass-init` are fixed in the script.

## Before re-running (see `TODO.md`)

- Both BayPass scripts set `-nval` (`NVAL=5000`; 3.1's default is 250).
  Expect roughly 25–50 h per job with 8 CPUs (2.31 runs with ≤ 35,000
  iterations took 9–12 h; 5,000 samples under 3.1 defaults is ~102,000
  iterations). That's 4 Ω jobs plus 20 covariate jobs. CHTC's default run
  limit is 72 h; the upper estimate is under it but not by much, so check
  the first job's progress before submitting the rest.
  BayPass 3.1 scales covariates by default (2.31 runs used raw logN).
- Existing Ω came from 2.31 with 8 threads, affected by a bug fixed in 2.4.
- Local analysis of outputs: `_R/baypass-check-omegas.R`,
  `_R/baypass-analysis.R`, `_R/baypass-check-candidates.R`.
