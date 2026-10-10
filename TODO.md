# TODO

Items 1–6 come from the 2026-10-09 code review (three independent reviewers;
key claims checked in `_R/review-checks.R` → `_data/review-checks.txt`).
On hold until the owner says to run them. Listed in dependency order.

## Before rerunning anything on the cluster

- [ ] **Rebuild the container** (`_bash/_apptainer/_build.sh`, v2.0.0). Pins
      were updated 2026-10-09 and dry-run solved for linux-64, but the image
      has not been built or tested. Outputs from the new image will differ
      from existing ones (old pins: `git show ca25cfc:_bash/_apptainer/`).
- [ ] **Re-estimate Ω.** Current Ω came from BayPass 2.31 with
      `-nthreads 8`; BayPass 2.4 fixed a bug in the Ω update "introduced with
      parallelisation in version 2.3x" (CHANGELOG). Effect size unknown, but
      all current candidates rest on this Ω. Rerun `01-baypass-init` with 3.1
      before any new covariate runs.
- [x] `-nval` set explicitly (`NVAL=5000`) in `01-baypass-init.sh` and
      `02-baypass-beta.sh` (3.1 default is 250). Expect ~25–50 h per job
      (estimate from 2.31 run times; see the scripts); CHTC's default run
      limit is 72 h. Covariate scaling is the
      3.1 default; noted in `02-baypass-beta.sh`. (2026-10-09)
- [x] `samtools mpileup -Q 25` in `03-mpileup.sh`, matching the counts'
      base-quality threshold. (2026-10-09)
- [x] `_R/baypass-utils.R` replaced with the 3.1 `utils/` version (local
      `concatenate_res` edits re-applied; listed in its header). (2026-10-10)
- [x] **SNP set:** `03-combine-syncs` now combines partially masked syncs
      and keeps sites SNAPE calls polymorphic (P ≥ 0.9) in at least one
      pool (option B). Previously a site had to be polymorphic in every
      pool, which excluded sites that drifted to fixation in any pool.
      (2026-10-10)
- [x] Dropped the `03-mpileup` window coverage (overlapping windows
      undercounted; output unused). (2026-10-10)
- [ ] With the new SNP set, `freq` in the SNAPE files is still each pool's
      own minor allele, and some pools are now fixed at a site. Review
      step 1 (one fixed allele per SNP from counts) matters more.

- [ ] **Pipeline options from DESTv3** (compared 2026-10-09; undecided):
      - DESTv3's `MaskSYNC_snape_monomorphic_filter.py` reassigns all reads
        to one allele at sites SNAPE calls confidently monomorphic (removes
        error reads). Low value for BayPass (SNP set is 0.1–0.9 mean freq).
      - Competitive mapping against other genomes (DESTv3 "hologenome") as a
        diagnostic of non-target reads per pool, if useful for review step 2.
        Cheaper first look: `02-bwa.sh` already saves unmapped reads
        (`*_unmapped_{merged,unmerged}.bam`) and `bamtools stats` output, so
        per-pool mapping rates vs. logN need no remapping. Other chironomid
        assemblies for a combined reference are in
        `~/_data/midgenomes/assemblies/`.

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
       - covariate scaling (manual §2.3.3 strongly recommends it; current
         2.31 runs used raw logN, mean 8.67, SD 1.59; scaled by default in
         3.1);
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

- [ ] **Pool inclusion** (decided 2026-10-10: QC-metric protocol instead of a
      median-coverage cutoff; implemented, not yet run).
      - Metrics: missing-data rate after masking, private SNPs, duplication
        rate, GC content, mapping rate (not pN/pS: annotation too poor).
        Computed by `_bash/02-allele-freqs/03-pool-qc`; bounds applied in
        `_R/pool-qc.R`. Set bounds after seeing the metrics but before any
        allele-frequency results, and record the rationale in the script.
      - Still to decide: fixed bounds per metric (DEST 1.0 style), or an
        outlier rule on the PCA of metrics (DEST 2.0 style).
      - Report the 2025 set (median coverage ≥ 30; `*-pools-cov30.txt`)
        as a sensitivity analysis.
      - Missing data per pool allowed in BayPass (masked pools → 0 counts;
        `00-baypass-inputs.sh`). **Provisional:** `MIN_PROP_POOLS = 0.75`
        (SNP kept if ≥ 75% of pools have data); confirm or change.
      - Background (`_R/coverage-cutoff-check.R` → `_data/coverage-cutoff-check.csv`):
        5 pools failed outright (median 0–4), including 1988a (N = 76),
        2009b, and 2013b; 1982a, 1988a, 1991b, and 2013b have GC 46–63%
        vs. 32.5% for the assembly. 8 pools have median 6–24, including
        1998b (N = 450). Conventions: DEST 1.0 (Kapun et al. 2021) kept
        pools down to 10x and excluded > 7% missing data; DEST 2.0 (Nunez et
        al. 2025) used a PCA of QC metrics; Schlötterer et al. 2014
        recommend ≥ 50x for SNP-level tests as a design target.
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
