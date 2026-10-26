# TODO

Items 1–6 come from the 2026-10-09 code review (three independent reviewers;
key claims checked in `_R/review-checks.R` → `_data/review-checks.txt`).
On hold until the owner says to run them. Listed in dependency order.

## Before rerunning anything on the cluster

Follow `_bash/RERUN-CHECKLIST.md`.

- [x] **Rebuild the container** (`_bash/_apptainer/_build.sh`) as v2.3.0:
      v2.2.0 (pins updated 2026-10-09) plus medaka model
      `r941_min_sup_g507` and `bin/mito-collapse.sh`; built by the owner
      2026-10-26 (RERUN-CHECKLIST C3). Outputs from the new image will differ
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
      `03-baypass-beta.sh` (3.1 default is 250). Expect ~25–50 h per job
      (estimate from 2.31 run times; see the scripts); CHTC's default run
      limit is 72 h. Covariate scaling is the
      3.1 default; noted in `03-baypass-beta.sh`. (2026-10-09)
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
        `00-references/02-hologenome` (T. gracilentus + *C. riparius*,
        *C. tentans* 2023, *Smittia aterrima*, *Propsilocerus akamusi*;
        decided 2026-10-12); `02-bwa` maps
        to it and counts reads per species; `03-mpileup` keeps only
        T. gracilentus contigs. Other species could be added (`OTHERS` in
        `00-references/02-hologenome.sh`); copy assemblies to `/staging/l/lnell/other-genomes/`
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
        contamination index. `00-references/01-mitogenome` assembles
        the T. gracilentus mitogenome from the nanopore reads of one
        identified male (R9.4.1, Dorado 0.9.6 sup v3.6 since 2026-10-23) with
        Flye + medaka (`r941_min_sup_g507`) + polypolish (pool SN-1-08_S24);
        check `coi-check.tsv` against COI MZ659759/MN679642 before going on.
        `00-references/02-hologenome` adds the mitogenome and COI barcodes
        (`_R/get-mito-refs.R` → `_bash/references/`; selection in
        `_data/coi-barcode-selection.csv`; *Procladius* at genus level since
        *P. islandicus* has no NCBI sequences). Index per pool and taxon:
        `_R/mito-contam.R` (relative, not a DNA fraction).
        - [ ] **`01-mitogenome` fails (first run 2026-10-14); decide how to
          build the reference.** Round-1 Flye assembles nothing at ~2,000x
          coverage (`--genome-size 16k --asm-coverage 150` fixes that).
          Every successful assembly is a 28-32 kb linear contig: ~13.6 kb
          collinear with Tanytarsus gene order plus AT-rich sequence with a
          ~2.8-kb tandem repeat where the control region sits; trimming,
          stricter recruitment, and `--meta` all fail `pick_contig`
          (13-22 kb). Long reads (> 17 kb) mostly match the gene region at
          mitochondrial identity, but 8 of 60 carry a short internal piece
          flanked by > 2 kb of other sequence (NUMT-like). Illumina DNA
          depth (pool MyKS-19-B, 4M pairs): the control-region-side
          sequence is at mitochondrial-level depth (140-2,400x vs 500-800x
          in the gene region; nuclear ~13x). Best supported: a long,
          repeat-rich control region (mitogenome >= ~22 kb); the circle is
          not closed and the repeat array's length is unreliable. The 2019
          NextDenovo genome assembly has a 66.4-kb contig with the gene
          region twice, 44,625 bp apart, separated by ~28.7 kb of control
          region including ~8 copies of a ~2.5-kb unit (copies ~89%
          identical); SMARTdenovo shows ~5 copies on each side of the gene
          region. If right, the mitogenome is ~44.6 kb (unverified).
          Goal (owner, 2026-10-14): a publication-quality mitogenome.
          - [x] Re-basecall all nanopore reads (`_bash/00-basecall/`; set up
            2026-10-14): archive on staging, `00-setup` and `01-split-fast5`
            run, `02-dorado` submitted 2026-10-22, `03-combine` done; local
            copy in `/Volumes/lan_images4/ont/` passed
            `_bash/tests/dorado-sup-check.sh` (2026-10-23: 3,602,836 reads
            at q >= 10, 23.0 Gb, N50 7.5 kb; median predicted q +3.1 over
            Guppy on the same reads). FAST5 archive and chunks deleted
            from staging by the owner after their checks (2026-10-26); the
            raw signal remains on `/Volumes/TanyGen5TB/nanopore/`.
            Then redo the mitochondrial recruitment and assemblies on the sup
            reads (owner chose to wait for them, 2026-10-19), and switch
            `01-mitogenome` to the sup FASTQ and medaka to a sup model
            (container rebuild to bundle it). `04-ont-pi` switched
            2026-10-23 (outputs tagged `_dorado-0.9.6_sup`; not yet run).
            Medaka has no model for `dna_r9.4.1_e8_sup@v3.6`; its
            basecaller table maps `dna_r9.4.1_e8_sup@v3.3` to
            `r941_min_sup_g507` (`medaka/options.py`, checked 2026-10-23).
            Owner to decide: re-basecall with v3.3, use `r941_min_sup_g507`
            on v3.6 reads, or polish with Illumina only. Local mitogenome
            tests rerun on sup reads from 2026-10-23 (Guppy cache kept in
            `_data/mito-flye-cache-guppy/`). Round 1 on sup reads
            (`mito-flye-round1.sh`): 7,221 reads recruited; default Flye
            still assembles nothing; with `--genome-size 16k --asm-coverage
            150`, one 24.2-kb linear contig (555x) that matches Guppy contig
            1,639-25,875 at 99.6% (BLAST), so `pick_contig` still fails.
            The mito checks now locate features by alignment
            (`_bash/tests/mito-features.sh`, 2026-10-23;
            `mito-features-check.sh` reproduces the old Guppy coordinates).
            Medaka model test (`medaka-model-test.sh` →
            `_data/medaka-model-test.txt`, sup round-1 contig vs MyKS-19-B
            Illumina): in the gene region, `r941_min_sup_g507` on the v3.6
            reads fixes all 6 of the draft's indel errors and adds 1;
            `r941_min_hac_g507` ends with 3; SNVs (17, likely
            pool-vs-individual differences) unchanged. Outside the gene
            region both models add the same 3 SNVs (sup one more, alt 0.91)
            where Illumina is not a reliable truth (repeat copies collapsed
            on the 24-kb contig; one site at 4,500x). Owner chose
            `r941_min_sup_g507` on the v3.6 reads (2026-10-23):
            `01-mitogenome` now reads the sup FASTQ and uses it; container
            v2.3.0 bundles it (rebuild needed, RERUN-CHECKLIST C3).
            Reassembly on sup reads (`mito-reassemble.sh`, 2026-10-23;
            2019 assemblies now read from `/Volumes/lan_images4/assembly/`):
            61,245 reads recruited. Flye's 26.1-kb contig and Raven's Utg158
            are nuclear (CONTIG12 63.8-89.2 kb at 99.7%, BLAST; not in the
            script output); Raven Utg156 (42.5 kb) has a gene region with
            two gaps and 26 kb of unidentified sequence; SMARTdenovo 53.6 kb
            with the gene region once and 9 unit copies, every boundary
            spanned by 25-140 reads; Illumina unit copies sum 4.2;
            NextDenovo failed at `seed_cns` again. `mito-cr-span.sh`: still
            no read crosses the control region (max 22.2-22.3 kb past
            either end). `mito-long-reads.sh`: same conclusions, identity
            ~0.5 points higher. `mito-depth-by-length.sh`: regions now
            differ (unit from another copy; 1 left, 8 right copies), so
            depths are not comparable with the Guppy run.
            Nuclear-read filter (2026-10-26): of the 61,245 recruited
            reads, ~3,600 align best to nuclear contigs (many small ones,
            CONTIG36-45; 62 long reads to CONTIG12), and among reads
            >= 11 kb, which Filtlong keeps, nuclear outnumber mitochondrial
            (1,050 vs 753). The assembly's own mitogenome is CONTIG34 (gene
            region twice; carries COI MZ659759.1). Added competitive
            mapping (`keep_mito`) to `01-mitogenome` recruitment and to
            `mito-reassemble.sh`. Rerun (2026-10-26): kept 57,864 reads,
            dropped 3,361 (1,117 of the 1,999 reads >= 11 kb); 20 reads with
            no alignment in the competitive mapping are dropped and not
            counted in `keep_mito.txt`. Flye now gives two mitochondrial
            contigs (26.5 and 28.4 kb, full gene region at >= 99.9%, 2-3
            unit copies); Raven gives no contigs; SMARTdenovo 34.2 kb with
            only 7.5 kb of the gene region (97.4%) and 7 unit copies (was
            53.6 kb, full gene region, 9 copies, unfiltered); NextDenovo
            fails at `seed_cns`. Subsample is now 1,966 reads down to 9.0 kb
            (was 1,480 down to 11.3 kb). The filter fixes Flye's nuclear
            contig but no assembly resolves the repeat array; no read spans
            the control region, so assembled lengths reflect each
            assembler's repeat handling.
          - Control-region length from reads (2026-10-26,
            `_bash/tests/mito-copy-number.sh` → `_data/mito-copy-number.txt`):
            no crossing reads give only a lower bound (5%: 18.2 kb; weaker
            than the 22.3-kb single-read extension). Nanopore depth can't
            count copies: along the 2019 NextDenovo contig (gene region +
            full control region), gene-region depth falls from ~22,000x
            (first 10 kb, rRNA side) to ~3,500x (ND2 side), and the control
            region is 208-640x except 8,385x in the 2.1 kb next to the
            gene region's start. The one-unit reference gives < 1 copy
            (0.3-0.7), contradicting reads spanning 7-8 unit boundaries.
            Cause of the uneven depth unknown (length heteroplasmy, partial
            molecules, A+T-rich read breakage, mapping loss).
          - Decided (owner, 2026-10-26): report the mitogenome as the gene
            region plus a control region of > 22 kb whose length can't be
            resolved with these data. `01-mitogenome` changed to match
            (2026-10-26): Flye with `--genome-size 16k --asm-coverage 150`
            in both rounds, the contig carrying COI MZ659759.1 is kept (no
            length rule), and after polishing `mito-collapse.sh`
            (`_apptainer/bin/`) keeps the gene region, the unique
            control-region sequence on each side, and one repeat unit
            (owner chose to include one unit). On the contigs at hand
            (`_bash/tests/mito-collapse-check.sh`) this gives 22.7-25.9 kb:
            genes 15.9 kb, unique right 3.5 kb (6.3 kb on 2019 SMARTdenovo),
            unit 2.37-2.48 kb, unique left 1.13 kb (0.97 kb, truncated, on
            the reassembled Flye contigs). Not run on the cluster (needs
            container v2.3.0). COI check now requires only MZ659759.1 at
            >= 98% (`REQUIRED_COI`); MN679642.1 (93.9-94.2% on every contig)
            is reported, not required (owner's decision, 2026-10-26). A
            possible reason for the gap (owner): MZ659759.1 is from Finland,
            MN679642.1 from Canada, so the Canadian barcode may come from a
            divergent lineage.
            Staging (2026-10-19): ~350 GB of 1,000 used after the owner
            removed the trimmed FASTQs (archive included); the split adds
            ~340 GB. `01-mitogenome` (polishing pool SN-1-08_S24) and
            `01-variant-calling/02-bwa` read the trimmed reads, so restore
            them (or rerun `01-fastp`) before those steps.
          - Reassembly (`mito-reassemble.sh`, rerun 2026-10-19): Flye
            38.4 kb (gene region partial; 6 repeat copies), Raven 42.8 kb
            (4), SMARTdenovo 52.4 kb (8), 2019 NextDenovo 66.4 kb (8, gene
            region twice 44.6 kb apart); NextDenovo failed at `seed_cns`.
            Every unit boundary is spanned by 19-420 uniquely aligned reads
            in all assemblies (not discriminating). Illumina (all reads)
            total unit copies ~4.6-5.6 in every assembly (NextDenovo's 7.4
            halves gene depth by having the gene region twice; ~4.7
            corrected); likely low, since Illumina under-covers A+T-rich
            sequence. With the > 22.5-kb control-region bound, 4 copies is
            too few; 6-8 copies (mitogenome ~42-45 kb) fit. Revisit with sup
            basecalls.
          - Control-region length (2026-10-14, `mito-cr-span.sh`): of 61,619
            mitochondrial reads, none crosses from the gene region's end to
            its start (wrap detection checked on synthetic reads); reads run
            up to 22.3-22.5 kb past either end without wrapping, so sampled
            molecules have a control region > ~22.5 kb (mitogenome
            >= ~38 kb). Nanopore depth falls with distance from the gene
            region (`mito-depth-by-length.sh`: 83% of gene-region depth in
            the 1.1 kb next to its start, 22% in the 3.6 kb next to its end);
            cause unresolved (repeat MAPQ, length variation among molecules,
            or assembly errors).
          (TanyAdult_S1/TanyJuven_S2 are RNA-seq of Myvatn adults and
          juveniles, so their ~0 depth there says nothing about it.)
          **COI check:** the DNA pool and the adult RNA-seq match MZ659759.1
          (98.6-98.8%) but MN679642.1 only at ~94.5% (4 and 767 reads vs
          2,665 and 220,905), so requiring >= 98% for both barcodes would
          reject a correct assembly. Scripts and outputs:
          `_bash/tests/mito-*.sh` → `_data/mito-*.txt`.
      - [x] Kraken2 diagnostic of unmapped reads (`01-variant-calling/05-kraken`;
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
      rRNA Kraken2 database (`00-references/03-silva-db`) run on
      all unmapped reads per pool in `01-variant-calling/05-kraken`, alongside PlusPF-16 on
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
