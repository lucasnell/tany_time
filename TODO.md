# TODO

Items 1–6 come from the 2026-10-09 code review (three independent reviewers;
key claims checked in `_R/review-checks.R` → `_data/review-checks.txt`).
On hold until the owner says to run them. Listed in dependency order.

## Review steps

1. [ ] **Rebuild `timeseries-data.csv` from counts** (`tany.genobaypass` +
       `tany.snpdet`) with one fixed allele per SNP and an allele column.
       Currently `freq` is the BayPass reference allele at 17 candidates and
       the alternative at 16. Small; local.
2. [ ] **Explain the shared direction of the candidates.** All 33 have the
       on-average-rarer allele more common at low abundance (random SNPs: 52%).
       Species contamination is ruled out (genitalia-checked males; no other
       *Tanytarsus* in Mývatn). Remaining explanations to test:
       - noise in the smallest, lowest-coverage pools (1997a/b, 2013a: 48, 38,
         28 haplotypes) combined with picking the tail of 433k tests;
       - a cryptic lineage or immigrants at low abundance;
       - few-family sampling in low-abundance years.
       Diagnostics vs. logN and pool size: mapping rate, coverage,
       heterozygosity in size-matched subsamples, PCA of pools. ~0.5 day;
       local.
3. [ ] **Leave-out BayPass runs:** without the 1997 pools (candidates' median
       |r| with logN drops 0.58 → 0.27 without them), and with generation
       (a/b) as a second covariate. Cheap; cluster.
4. [ ] **Calibration and chain length, reported as a full grid:**
       - pseudo-observed data (POD) via `simulate.baypass()` (manual §4.1)
         with real Ω, pool sizes, and coverages, to set a null BF threshold;
       - longer chains (`-nval`); single runs currently give 178–235 SNPs
         ≥ 30 dB and 1,304 distinct across 10 runs; 52.96 dB is a PIP = 1
         ceiling at 1,000 samples;
       - `-scalecov` (manual §2.3.3 strongly recommends scaling; current
         runs use raw logN, mean 8.67, SD 1.59);
       - Ising β = 0 and 1 (β = 1 only links SNPs adjacent in the input file;
         ~44 contig junctions);
       - covariates: logN, lagged logN, and a crash-specific covariate
         (crashes are the main question).
       Report all settings, not the best one. 1–2 days; cluster.
5. [ ] **Neutral variation through time:** read the `02-npstat` output
       (currently unused in `_R/`); per-pool diversity with pool-size-aware
       estimators; Ne through time from allele-frequency change between
       consecutive pools. Use an unascertained SNP set (the BayPass set drops
       mean frequency outside 0.1–0.9). Independent of 1–4. 3–5 days; local.
6. [ ] **Drift vs. selection:** test whether candidate frequency changes
       exceed what genome-wide Ne implies. After 2–4.

## Other analysis items

- [ ] **Pool coverage cutoff.** Median coverage ≥ 30 drops 13 of 30 SN pools,
      including 1988a (N = 76), the only sequenced deep-crash pool. Reason
      for the cutoff is not recorded; owner is open to changing it. Revisit
      with reports in `~/_data/_tany-time/trimmed/trimmed_reports/`. Since
      no new pools can be added, this is the only way to widen the abundance
      range.
- [ ] GO enrichment uses the default universe (all annotated *D.
      melanogaster* proteins) with 8 input genes; unusable until the
      universe is genes near tested SNPs with a Dmel ortholog. Count linked
      candidates as one locus first (CONTIG43 308275/308276 are 4 bp apart;
      CONTIG23 36355/36363 ~1 kb).
- [ ] No upper-coverage SNP filter (paralogs/collapsed repeats can cause
      frequency swings); `popsync2pooldata` uses default coverage settings.
- [ ] Check `n_eff` formula in `make-time-series-data.R` against Feder et al.
      2012 (one reviewer questioned the denominator).
- [ ] Time-series analysis for balancing selection — after 1–6.

## Housekeeping

- [ ] `_R/sample-qc.R` points to `~/_data/trimmed_reports`; reports are now in
      `~/_data/_tany-time/trimmed/trimmed_reports/`.
- [ ] `_R/make-time-series-plot.R` saves its main figure to
      `~/Desktop/freq-abund-time.pdf`; point it at `_plots/`.
- [ ] Move BLAST inputs/outputs off `~/Desktop` (`Tgraci_proteins.faa`,
      `Tgraci-funct-Dmel.csv`) to a stable path.
- [ ] `baypass-check-omegas.R`: the comment `range(fmd_df$fmd)` → 0.30–0.38
      contradicts the printed table (max 3.78). Update the stale comment.
- [ ] `candi_df.rds`, `candi_funcs_df.rds`, and `timeseries-data.csv` are
      written only if absent; stale after any rerun.
- [ ] `README.md` contains only "midge_archive".
