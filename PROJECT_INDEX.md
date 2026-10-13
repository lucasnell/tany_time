# Project index

## Pipeline order

Cluster steps are documented in `_bash/README.md` and each `_bash/*/README.md`;
`_bash/RERUN-CHECKLIST.md` is the step-by-step list for a full rerun.
Job lists (`_bash/samples/`, `_bash/03-select-time/baypass_beta_info.txt`)
are written by `_R/make-job-lists.R`.

1. `_bash/01-variant-calling/` (cluster): `00a-mitogenome` (T. gracilentus
   mitogenome from nanopore reads) → `00b-hologenome` (combined reference:
   T. gracilentus + 2 *Chironomus* + 2 Orthocladiinae nuclear genomes, the
   mitogenome, and COI barcodes of trap taxa) → `01-fastp` trim →
   `02-bwa` merge overlapping pairs with bbmerge, align competitively with
   `bwa mem` → `03-mpileup` mark duplicates, GATK indel realignment,
   `samtools mpileup` on T. gracilentus contigs; `05-bam-cov` (coverage
   `.bed.gz`); `06a-silva-db` + `06b-kraken` (Kraken2 on unmapped reads: PlusPF and SILVA rRNA; diagnostic);
   `04-ont-pi` (SNAPE prior θ, independent). Competitive mapping and
   Kraken2 added 2026-10-10; existing outputs mapped to T. gracilentus only.
2. `_R/tidy-sample-data.R` → `_data/full-sample-info.csv` (pool ID, date,
   lake, site, number of adults, lat/lon).
3. `_R/sample-qc.R`: picks pools with median coverage ≥ 30 (17 of 30 SN time
   pools, 18 spatial pools; lists are in the script). Reason for the cutoff
   is not recorded. The script's report path (`~/_data/trimmed_reports`) is
   stale; reports are now in `~/_data/_tany-time/trimmed/trimmed_reports/`.
4. `_bash/02-allele-freqs/` (cluster; SNAPE helper scripts, now in
   `_bash/_apptainer/bin/`, come from DEST freeze1 and match DESTv3 except
   for uncompressed output and the P(polymorphic) filter in
   `MaskSYNC_snape.py`, which keeps only sites with P ≥ 0.9; DESTv3 masks
   only 0.1 < P < 0.9): `01-snape` (SNAPE-pooled → masked gSYNC; undoes
   softmasking) → `02-npstat`; `03-pool-qc` (per-pool QC metrics) →
   `_R/pool-qc.R` (locally; chooses pools) → `04-combine-syncs` (combines included pools
   into `snape_masked_*` and `count_masked_*` for `time`, `space`, and `all`;
   was `03-combine-snapes` + `04-combine-counts` before 2026-10-09). Since
   2026-10-10 it uses partially masked syncs and keeps sites polymorphic in
   at least one pool; existing outputs required polymorphism in every pool.
   Since 2026-10-13, count files assign all reads to the reference in pools
   SNAPE calls confidently monomorphic (P ≤ 0.1).
   BayPass inputs now allow missing data per pool (masked pools → 0 counts).
5. `_R/tidy-pop-abunds.R` → `_data/tany-abundances.csv` (year, generation,
   N, logN = log1p(N)).
6. `_bash/03-select-time/` (cluster):
   - `00-baypass-inputs`: poolfstat filtering, drops SNPs with mean
     reference-allele read proportion outside (0.1, 0.9), writes genotype,
     pool-size, and logN covariate files. 433,360 SNPs.
   - `01-baypass-init`: Ω (population covariance matrix) from all SNPs
     (`sub0`) and 3 subsamples. Existing Ω is from BayPass 2.31 with 8
     threads, affected by an Ω-update bug fixed in 2.4 (see `TODO.md`).
   - All three steps run for a grid of missing-data settings
     (`MIN_PROP_POOLS` 0.5, 0.75, 1.0; tags `mp050`/`mp075`/`mp100`;
     0.75 main; since 2026-10-13).
   - `02-baypass-beta`: covariate model with Ω from all SNPs; 10 runs with
     different seeds; Ising β given as an argument. Current runs use
     `-auxmodel` with β = 0 (no dependency between neighbouring SNPs),
     because the assembly is not chromosome-level (contig-to-contig
     linkage is unknown). 1,000 retained MCMC samples per run. The covariate
     is raw logN (not `-scalecov`), though the BayPass manual (§2.3.3)
     strongly recommends scaling.
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
| `baypass-utils.R` | BayPass v3.1 R utilities (`simulate.baypass`, `fmd.dist`, `concatenate_res`, …) with local `concatenate_res` edits |
| `baypass-check-omegas.R` | Ω consistency across subsamples |
| `baypass-analysis.R` | Candidate SNPs, gene matching, BLAST, enrichment |
| `baypass-check-candidates.R` | Confirms `candi_df.rds` matches local BayPass output; per-candidate run consistency → `_data/baypass-check-candidates.txt` |
| `make-time-series-data.R` | Builds `timeseries-data.csv` |
| `make-time-series-plot.R` | Candidate frequency vs. abundance plots |
| `make-job-lists.R` | Writes HTCondor job lists in `_bash/samples/` (all pools, candidate pools, coverage-30 sets) and `_bash/03-select-time/baypass_beta_info.txt` (reproduces the 2025 BayPass seeds) |
| `pool-qc.R` | Chooses pools from `03-pool-qc` metrics: ROBPCA outlier rule with bootstrap-calibrated threshold (blind to abundance and allele frequencies) → `_bash/samples/{time,space}-pools.txt`, `_data/pool-qc.csv`, `_plots/pool-qc-*.pdf`; needs `rrcov` |
| `pool-qc-rule-sim.R` | Simulates false-positive and detection rates of outlier-rule variants on clean and planted data → `_data/pool-qc-rule-sim.csv`, `_data/pool-qc-rule-sim-calibrated.csv` |
| `coverage-cutoff-check.R` | Per-pool coverage, duplication, and GC for SN pools with abundance → `_data/coverage-cutoff-check.csv` |
| `get-mito-refs.R` | Fetches T. gracilentus mitochondrial bait/COI and trap-taxon COI barcodes from NCBI → `_bash/references/*.fasta`, `_data/coi-barcode-selection.csv` |
| `kraken-summary.R` | Summarizes Kraken2 reports of unmapped reads (`06b-kraken`): PlusPF composition (bacteria, fungi, human, etc.; top genera; plot) and SILVA parasite groups relative to host rRNA (top non-insect eukaryotes) → `_data/kraken-groups.csv`, `_data/kraken-top-genera.csv`, `_plots/kraken-groups.pdf`, `_data/kraken-silva-groups.csv`, `_data/kraken-silva-top-eukaryotes.csv`; not yet run |
| `mito-contam.R` | Per-pool mitochondrial contamination index per trap taxon from `02-bwa` read counts → `_data/mito-contam.csv`; not yet run |
| `review-checks.R` | Checks of the 2026-10 code-review claims (polarity, shared direction, 1997 dependence, per-run counts, crash sampling) → `_data/review-checks.txt` |

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

From `_R/review-checks.R` → `_data/review-checks.txt` (run 2026-10-09):

- `freq` in `timeseries-data.csv` tracks the BayPass reference allele at 17
  candidates and the alternative allele at 16; it is not a consistent allele.
- Coded by the allele rarer on average across pools, 33 of 33 candidates have
  that allele more common at low abundance; 52.3% of 100,000 random SNPs do
  (binomial P = 5e-10).
- Without the two 1997 pools (lowest and third-lowest logN), median
  |Pearson r| of candidates with logN drops from 0.58 to 0.27; 6 of 33 stay
  nominally significant (n = 15).
- Single runs give 178–235 SNPs ≥ 30 dB; 1,304 distinct SNPs pass in at
  least one run.
- 0.35% of random SNPs have |r| > 0.7 with logN vs. 0.18% expected for
  independent tests.
- 13 of 30 sequenced SN pools fail the coverage cutoff, including 1988a
  (N = 76). Minimum N among used pools is 246; 17 generations in 1977–2015
  had N < 100.

From `_data/candi_funcs_df.rds`: 26 of 33 candidates fall inside or within
5 kb of a gene (24 unique genes); 8 unique *D. melanogaster* UniProt best hits
at e < 1e-4 (listed with annotations in `baypass-analysis.R`).

## `_bash/tests/`

Local checks of pipeline behavior: `npstat-maxcov.sh` (lowering npstat's
`-maxcov` to the observed maximum leaves output unchanged; →
`_data/npstat-maxcov-check.txt`).

## `_refs/`

`BayPass_manual_v2.31.pdf` (version used for existing results) and
`BayPass_manual_v3.1.pdf` (version in the v2.0.0 container). Local copies;
PDFs are git-ignored.

## `_legacy/`

Older analyses moved from the repo root (2026-10-09); not maintained and not
runnable as-is (see `CLAUDE.md`). Includes sample summaries, spatial FST and
connectivity (`fst-dist.stan`), admixture, temporal diversity, diversity
simulations (`sim_midges.cpp` with the `pcg/` RNG headers), and genome-report
maps. Git-ignored outputs `fst-dist.rds` and `sample_summary.pdf` moved with
them. `cpg-islands.sh` (local TaJoCGI run on the assembly; unfinished, found
only 7 CpG islands) moved here from `_bash/` on 2026-10-09.
