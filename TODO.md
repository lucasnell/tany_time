# TODO

Items 1–6 come from the 2026-10-09 code review (three independent reviewers;
key claims checked in `_R/review-checks.R` → `_data/review-checks.txt`).
On hold until the owner says to run them. Listed in dependency order.

## Before rerunning anything on the cluster

- [ ] **Rebuild the container** (`_bash/_apptainer/_build.sh`, v2.2.0). Pins
      were updated 2026-10-09 and dry-run solved for linux-64, but the image
      has not been built or tested. Outputs from the new image will differ
      from existing ones (old pins: `git show ca25cfc:_bash/_apptainer/`).
- [ ] Update `_R/baypass-analysis.R`, `baypass-check-omegas.R`, and
      `baypass-check-candidates.R` for the new BayPass output names
      (missing-data tag) and to compare settings.
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

- [x] **Pipeline options from DESTv3** (compared 2026-10-09; all decided):
      - DESTv3's `MaskSYNC_snape_monomorphic_filter.py` reassigns all reads
        to one allele where SNAPE calls a pool confidently monomorphic
        (P ≤ 0.1) or fixed for the alternative (removes error reads). More
        relevant since option B (2026-10-10): BayPass sites can now include
        pools fixed at that site, whose raw counts contain error reads.
        Cost: erases true rare alleles that SNAPE calls monomorphic.
        Decided 2026-10-13: narrow version only (P ≤ 0.1 → all reads to
        the reference; no rewrite for pools fixed for the alternative),
        `--clean-mono` in `combine-syncs.py`, used for count files in
        `04-combine-syncs` (`MAX_PP_MONO=0.1`).
      - [x] Competitive mapping (DESTv3 "hologenome") added 2026-10-10:
        `01-variant-calling/00b-hologenome` (T. gracilentus + *C. riparius*,
        *C. tentans* 2023, *Smittia aterrima*, *Propsilocerus akamusi*;
        decided 2026-10-12); `02-bwa` maps
        to it and counts reads per species; `03-mpileup` keeps only
        T. gracilentus contigs. Other species could be added (`OTHERS` in
        `00b-hologenome.sh`); copy assemblies to `/staging/lnell/other-genomes/`
        first.
        Species chosen from Einarsson et
        al. 2026 (Table 1) and NCBI (queried 2026-10-12): assemblies exist
        for *Chironomus* (incl. chromosome-level *C. tentans*
        GCA_963573255.1, better than the local 2014 scaffolds) and for
        Orthocladiinae (*Smittia aterrima* GCA_033063855.1, present but rare
        at Mývatn; *Propsilocerus akamusi*, local), but none for
        *Cricotopus*, *Micropsectra* (Tanytarsini, the closest relatives of
        T. gracilentus), *Procladius*, *Eukiefferiella*, or *Simulium*.
      - [ ] **Micropsectra (Tanytarsini) contamination** (discussing
        2026-10-12). NCBI (queried 2026-10-12) has no Micropsectra or other
        non-Tanytarsus Tanytarsini genome, transcriptome, raw reads, or
        mitogenome; it has ~4,260 Micropsectra COI records (including
        *M. lindrothi* voucher IS39, KC788724) and Tanytarsini mitogenomes
        (*Tanytarsus formosanus*, *T. verralli*, *Cladotanytarsus
        pseudomancus*). The T. gracilentus assembly appears to lack a
        mitogenome (CONTIG1, 16 kb, has no match to the T. gracilentus COI
        barcode MZ659759). Decided and implemented 2026-10-12 (not yet run):
        mitochondrial/barcode competitive mapping as a per-pool
        contamination index. `01-variant-calling/00a-mitogenome` assembles
        the T. gracilentus mitogenome from the nanopore reads of one
        identified male (R9.4.1, Guppy 5.0.11 hac; not re-basecalled) with
        Flye + medaka (`r941_min_hac_g507`) + polypolish (pool SN-1-08_S24);
        check `coi-check.tsv` against COI MZ659759/MN679642 before going on.
        `00b-hologenome` adds the mitogenome and COI barcodes
        (`_R/get-mito-refs.R` → `_bash/references/`; selection in
        `_data/coi-barcode-selection.csv`; *Procladius* at genus level since
        *P. islandicus* has no NCBI sequences). Index per pool and taxon:
        `_R/mito-contam.R` (relative, not a DNA fraction).
      - [x] Kraken2 diagnostic of unmapped reads (`01-variant-calling/06b-kraken`;
        `kraken-env` in container v2.1.0+). Needs a Kraken2 database on
        staging (PlusPF-16 suggested). Reports summarized by
        `_R/kraken-summary.R` (key groups, top genera, plot).
        Cheaper first look: `02-bwa.sh` already saves unmapped reads
        (`*_unmapped_{merged,unmerged}.bam`) and `bamtools stats` output, so
        per-pool mapping rates vs. logN need no remapping. Other chironomid
        assemblies for a combined reference are in
        `~/_data/midgenomes/assemblies/`.

- [x] **npstat speed** (2026-10-13): npstat precomputes tables for every
      coverage up to `-maxcov` at a cost ~ n × maxcov² (binomial terms), and
      `02-npstat` passed `-maxcov 10000` for each of ~45 contigs per pool.
      Now uses each contig's maximum coverage. Output identical, ~90x
      faster on a 200-kb test (`_bash/tests/npstat-maxcov.sh` →
      `_data/npstat-maxcov-check.txt`).
- [ ] **Parasites with Kraken2** (discussing 2026-10-13). Sources:
      `~/Downloads/midgeparasites/` (2015 Mývatn parasite protocols and
      literature, incl. de Souza et al. 2014 review of fungal/oomycete
      parasites of Chironomidae). Prebuilt Kraken2 databases are
      RefSeq-based: they include endosymbiotic bacteria (*Wolbachia*,
      *Spiroplasma*) and many viruses, but none of the chironomid parasite
      genomes, which are GenBank-only (`GCA`): *Entomophthora muscae*,
      *Coelomomyces lativittatus*, *Catenaria anguillulae*, *Lagenidium
      giganteum*, *Romanomermis*/*Mermis* (mermithids). No genomes for
      water mites (Hydrachnidia), *Aphanomycopsis*, *Atkinsiella*,
      *Couchia*, *Crypticola*, *Culicinomyces*, or *Amblyospora*
      (NCBI, queried 2026-10-13). Adults (our pools) can carry
      *Entomophthora*, water-mite larvae, mermithids, vertically
      transmitted microsporidia, endosymbiotic bacteria, and viruses;
      egg and larval parasites are unlikely in adults.
      Decided and implemented 2026-10-13 (not yet run): a SILVA 138.2 SSU
      rRNA Kraken2 database (`01-variant-calling/06a-silva-db`) run on
      all unmapped reads per pool in `06b-kraken`, alongside PlusPF-16 on
      1M reads; summarized per parasite group, relative to host
      (Chironomidae) rRNA reads, by `_R/kraken-summary.R`. Check its
      warning about SILVA names that matched nothing. Caveat: parasite DNA
      may come from other specimens in the trap preservative.

## Review steps

1. [ ] **Rebuild `timeseries-data.csv` from counts** (`tany.genobaypass` +
       `tany.snpdet`) with one fixed allele per SNP and an allele column.
       Currently `freq` is the BayPass reference allele at 17 candidates and
       the alternative at 16. Small; local.
2. [ ] **Explain the shared direction of the candidates.** All 33 have the
       on-average-rarer allele more common at low abundance (random SNPs: 52%).
       Misidentified *Tanytarsus* is ruled out (genitalia-checked males; no
       other *Tanytarsus* in Mývatn). DNA carried over from other trap taxa
       (more likely when T. gracilentus was scarce) is now tested by
       competitive mapping, the mitochondrial contamination index
       (`_R/mito-contam.R`), and Kraken2 (`_R/kraken-summary.R`). Other
       explanations to test:
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
        Computed by `_bash/02-allele-freqs/03-pool-qc`; rule applied in
        `_R/pool-qc.R` before any allele-frequency results exist.
      - Rule: outlier in a robust PCA of the metrics (decided 2026-10-10).
        Implemented as ROBPCA (`rrcov::PcaHubert`; `rrcov` 1.7.7 installed
        2026-10-10) with the outlier-score threshold calibrated by
        parametric bootstrap for a 2.5% per-pool false-positive rate.
        Uncalibrated versions, including ROBPCA's own cutoffs, excluded
        2.6–7 clean pools per set in simulations (`_R/pool-qc-rule-sim.R`
        → `_data/pool-qc-rule-sim*.csv`). The mapping-rate metric is
        `target_rate` (reads assigned to T. gracilentus under competitive
        mapping). Calibrated rule in simulation (200 reps;
        `_data/pool-qc-rule-sim-calibrated.csv`): false-positive rate 2.4%
        per pool; gross outliers (6 SD on 2 metrics) caught 94% (18 pools)
        and 99% (30 pools); a moderate outlier (4 SD on 1 metric) caught
        30% and 61%.
      - Report the 2025 set (median coverage ≥ 30; `*-pools-cov30.txt`)
        as a sensitivity analysis.
      - Missing data per pool allowed in BayPass (masked pools → 0 counts;
        `00-baypass-inputs.sh`); `MIN_PROP_POOLS` = minimum proportion of
        pools with data at a SNP. Literature search (2026-10-13): no standard
        value for pooled/population-level data; practice ranges from no
        missing data (e.g., Machado et al. 2021 required polymorphism in
        every population) to call-rate thresholds like ≥ 90% per
        population (kelp temporal study) and 0.6–0.8 in individual-level
        RAD/GEA work. Ahrens et al. (F-word paper) found GEA results are
        sensitive to small changes in missing-data filters; stricter
        filters cut power. Missingness here tracks coverage, so it isn't
        random with respect to abundance. Decided 2026-10-13: grid of 0.5,
        0.75, 1.0 (tags `mp050`, `mp075`, `mp100`), all reported, 0.75 main;
        implemented in `_bash/03-select-time/` (6 Ω jobs, 60 covariate
        jobs).
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

Done 2026-10-13: `sample-qc.R` report path; `make-time-series-plot.R` saves
to `_plots/`; BLAST files moved to `~/_data/_tany-time/blast/`;
`baypass-check-omegas.R` range comment (0.27–3.78); derived datasets
written unconditionally; repo `README.md` written.
