# Cluster rerun checklist

Full rerun of the pipeline with container v2.4.0 (v2.3.0 with bbmap back
to 39.28; v2.3.0 is v2.2.0 plus the sup medaka model). Paths are on CHTC
staging (`/staging/l/lnell/`, written `$S` below) unless noted. Submit each
step from inside its folder after `mkdir -p logs/<step>`. Details for each
step are in that folder's `README.md`.

Steps marked *(parallel)* can run alongside the step before them.

## A. Local preparation

- [x] **A1.** Pull the latest repo (`git pull`) locally and on the submit
      server (`/home/lnell/tany_time`); check `git status` is clean in both.
- [x] **A2.** Regenerate job lists (only needed if samples, QC, or the
      BayPass grid changed): `Rscript _R/make-job-lists.R`.
- [x] **A3.** Make the sample-info tar for BayPass and copy it to
      staging (commands at the top of
      `03-select-time/00-baypass-inputs.sh`) → `$S/tany-sample-info.tar.gz`.

## B. Staging setup (on the transfer/submit server)

- [x] **B1.** Create output folders (only pool QC, Kraken2, and the SILVA
      build create their own):

      ```bash
      S=/staging/l/lnell
      mkdir -p $S/dna/{raw_fq,trimmed,bwa,mpileup,ont,snape,npstat,pool-qc,kraken} \
          $S/dna/baypass/{init,beta} $S/{other-genomes,silva,kraken2-db}
      ```

- [x] **B2.** Confirm inputs already on staging:
  - [x] Raw reads: `$S/dna/raw_fq/<pool>_L002_R{1,2}_001.fastq.gz` for all
        62 pools in `samples/all-pools.txt`:

        ```bash
        while IFS= read -r line; do
          FILE=$S/dna/raw_fq/${line}_L002_R1_001.fastq.gz
          [ ! -f "$FILE" ] && echo "$FILE is missing."
          FILE=$S/dna/raw_fq/${line}_L002_R2_001.fastq.gz
          [ ! -f "$FILE" ] && echo "$FILE is missing."
        done < "/home/lnell/tany_time/_bash/samples/all-pools.txt"
        ```

  - [x] Nanopore reads: `$S/dna/ont/basecalls_guppy-5.0.11.fastq.gz`.
  - [x] Assembly: `$S/Tgraci_assembly.fasta.gz`.
  - [x] Repeats: `$S/Tgraci_repeats_locs.gff3.gz`.
- [x] **B3.** Copy the other assemblies from
      `~/_data/midgenomes/assemblies/` to `$S/other-genomes/`:
      `Cripar_assembly.fasta.gz`, `Ctenta2023_assembly.fasta.gz`,
      `Saterr_assembly.fasta.gz`, `Pakamu_assembly.fasta.gz`:

      ```bash
      cd ~/_data/midgenomes/assemblies
      scp Cripar_assembly.fasta.gz Ctenta2023_assembly.fasta.gz \
        Saterr_assembly.fasta.gz Pakamu_assembly.fasta.gz \
        transfer.chtc:/staging/l/lnell/other-genomes/
      ```

- [x] **B4.** Download the three SILVA 138.2 files into `$S/silva/` (URLs at
      the top of `00-references/03-silva-db.sh`):

      ```bash
      # on cluster:
      cd /staging/l/lnell/silva/
      wget https://ftp.arb-silva.de/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
      wget https://ftp.arb-silva.de/release_138_2/Exports/taxonomy/tax_slv_ssu_138.2.acc_taxid.gz
      wget https://ftp.arb-silva.de/release_138_2/Exports/taxonomy/tax_slv_ssu_138.2.txt.gz
      ```

- [x] **B5.** Download Kraken2 PlusPF-16 into `$S/kraken2-db/` from
      https://benlangmead.github.io/aws-indexes/k2, and set
      `KRAKEN_DB_TAR` in `01-variant-calling/05-kraken.sh` to its file
      name (currently a placeholder, `k2_pluspf_16gb.tar.gz`):

      ```bash
      # on cluster:
      cd /staging/l/lnell/kraken2-db/
      wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16_GB_20260626.tar.gz
      ```


- [x] **B6.** Clone the repo on the submit server (once) and pull before
      each session: `cd /home/lnell && git clone
      https://github.com/lucasnell/tany_time.git`, then
      `cd /home/lnell/tany_time && git pull`. Submit jobs from inside the
      clone; `logs/` folders are git-ignored.

## C. Container

- [x] **C1.** Build the image (`_apptainer/README.md`): copy
      `_apptainer/` files to the submit server, run `_build.sh` in an
      interactive build job → `$S/tany_time.sif`.
- [x] **C2.** Quick check in an interactive job: each environment
      activates (`main-env`, `qc-env`, `gatk-env`, `gatk3-env`,
      `baypass-env`, `kraken-env`, `mito-env`), `/app/bin` scripts are on
      `PATH`, and `medaka tools list_models` (in `mito-env`) shows
      `r941_min_hac_g507` (v2.2.0):

      ```bash
      bash
      . /app/.bashrc
      for f in $(find /app/bin -name "*.py"); do
        which $(basename $f)
      done
      conda activate main-env
      conda deactivate
      conda activate qc-env
      conda deactivate
      conda activate gatk-env
      conda deactivate
      conda activate gatk3-env
      conda deactivate
      conda activate baypass-env
      conda deactivate
      conda activate kraken-env
      conda deactivate
      conda activate mito-env
      medaka tools list_models | grep  --color=always r941_min_hac_g507
      ```

- [x] **C3.** Rebuild as v2.3.0 (medaka model `r941_min_sup_g507` for the
      sup reads, 2026-10-23; `bin/mito-collapse.sh`, 2026-10-26) as in C1,
      then C2 with `grep r941_min_sup_g507` and `which mito-collapse.sh`.
      Needed before F1.
- [x] **C4.** Rebuild as v2.4.0 (bbmap 39.28: 40.02's BBMerge dropped or
      hung on paired reads, so read mapping (G1) produced near-empty BAMs
      or stalled; 2026-10-29) as in C1, then rerun the tool check or
      run `bbmerge.sh --version` (expect 39.28). Needed before G1.

## D. How to run and check each step

From the step's folder on the submit server: `mkdir -p logs/<step>`, then
`condor_submit <step>.sub`. When the jobs finish (or to check on them while
they run), run `check-step.sh` from the `_bash/` folder:

```bash
cd /home/lnell/tany_time/_bash
./check-step.sh <folder>/<step>
```

It reports finished, failed, running, and unsubmitted jobs; the longest
run time (warning above 54 h; CHTC's limit is 72 h); outputs missing from
staging; and a step-specific summary. Go on only after `RESULT: PASS`.
For steps marked "one first", submit a single job and check its run time
and disk use before submitting the rest. Pick the largest input (worst
case). `condor_submit -queue` fails on a submit file that has its own
`queue` line, so rewrite that line and pipe the result in, e.g. for
`03-mpileup` (the `.sub` file is unchanged):

```
sed 's/^queue readbase from .*/queue readbase in (Lei-19_S13)/' 03-mpileup.sub | condor_submit -
grep -v '^Lei-19_S13$' ../samples/all-pools.txt > logs/03-mpileup/rest.txt
sed 's|^queue readbase from .*|queue readbase from logs/03-mpileup/rest.txt|' 03-mpileup.sub | condor_submit -
```

## E. Trimming (`01-variant-calling/`)

- [x] **E1.** `01-fastp` (62 jobs) → `./check-step.sh 01-variant-calling/01-fastp`.
      The trimmed reads (`.tar`) were deleted from staging on 2026-10-30
      after G1 (copies on lan_images4), so this check now reports them
      missing.

## F. References (`00-references/`)

- [x] **F1.** `01-mitogenome` (1 job; needs trimmed reads of
      SN-1-08_S24 from E1; PASS 2026-10-28 with the `mito-collapse.sh`
      fix and round-2 `--genome-size 45k`). It fails by itself, writing
      nothing, unless
      COI barcode MZ659759.1 matches the collapsed mitogenome
      (gene region, unique control-region sides, one repeat unit) →
      `./check-step.sh 00-references/01-mitogenome`.
- [x] **F2.** `02-hologenome` (1 job; drops CONTIG34, the assembly's own
      mitogenome, listed in the `01-mitogenome` output; other species'
      mitogenomes as `mt_<code>__`. PASS 2026-10-28, 0.2 h: 44 target
      contigs, 91,778,304 bp; 17 groups). It fails by itself unless every
      species, the mitogenome, and every barcode group are present →
      `./check-step.sh 00-references/02-hologenome`.
- [x] **F3.** *(parallel with F1–F2)* `03-silva-db` (1 job) →
      `./check-step.sh 00-references/03-silva-db`.
- [x] **F4.** *(optional, parallel)* `04-ont-pi` (1 job) →
      `./check-step.sh 00-references/04-ont-pi` compares the new θ with
      the SNAPE prior in `02-allele-freqs/01-snape.sh` and warns if they
      differ by more than 10%. Rerun on the sup reads 2026-10-23:
      0.008894 (+1.1% vs the old prior 0.0088); THETA in `01-snape.sh` set
      to this value, unrounded, 2026-10-28.

## G. Reads to pileups (`01-variant-calling/`)

- [x] **G1.** `02-bwa` (62 jobs) → `./check-step.sh 01-variant-calling/02-bwa`
      (also tabulates each pool's reads on *T. gracilentus*, other
      species, the mitogenome, and barcodes). First run (2026-10-28,
      container v2.3.0) failed: bbmap 40.02's BBMerge kept ~1-2% of read
      pairs and exited 0 (38 jobs; near-empty BAMs on staging) or hung
      with no CPU use (24 jobs). Remove those jobs, rebuild (C4), and
      rerun all 62; the job now fails if BBMerge loses any pair. Rerun
      submitted 2026-10-29 (v2.4.0, bbmap 39.28): all 62 finished OK
      (2026-10-30; longest 3.7 h), none held, every pool kept all its
      read pairs, 124 of 124 outputs on staging. Per-pool reads by group:
      `_data/02-bwa-check.txt`. The *T. gracilentus* mitogenome takes
      over half of the mapped reads in 6 pools (SN-1-96 88%, SN-2-96
      85%, SN-2-98 75%, H-2-06 72%, KS-1-15-B 57%, SN-1-13 56%); other
      species' nuclear genomes take 43% in SN-2-09 and 20% in H-2-10.
- [x] **G2.** `03-mpileup` (62 jobs, 9+ h each; one first) →
      `./check-step.sh 01-variant-calling/03-mpileup` (also tabulates
      reads after duplicate removal and after realignment, which must
      match, and pileup sizes; added 2026-10-31). Test job
      (Lei-19_S13, largest BAM at 9.6 GB) finished 2026-10-30: 4.7 h,
      3.8 CPUs, 9.0 of 32 GB memory, 45.8 of 50 GB disk; duplicates 19%
      (109.7 million records after removal, same after realignment); no
      errors. Outputs 10.2 GB (pileup 3.3 GB, tar.gz with the realigned
      BAM 6.8 GB). Scaled by BAM size, all 62 need ~280 GB; staging held
      719 of 1,000 GB; the owner deleted `$S/dna/trimmed/` (359 GB;
      only `02-bwa` reads the trimmed reads) and the other 61 were
      submitted 2026-10-30 (cluster 12392440) with `request_disk` raised
      to 60 GB. The deletion also removed the `trimmed_<pool>.tar.gz`
      fastp reports, which `03-pool-qc` (H3) reads for `fastp.json`
      (GC content is a rule metric); the owner restored them from
      lan_images4 the same day (62 files, all fastp 1.4.0, readable).
      All 62 finished OK (2026-10-31): PASS, longest 4.7 h (Lei-19_S13),
      124 of 124 outputs on staging (checked before the read-count
      table was added).
- [x] **G3.** *(parallel with G2)* `04-bam-cov` (62 jobs) →
      `./check-step.sh 01-variant-calling/04-bam-cov`. Run 2026-10-30
      (cluster 12391034): PASS, longest 0.4 h. Spot check of the smallest
      and largest outputs (KS-2-83, Lei-19): each covers all 44
      *T. gracilentus* contigs end to end (91,778,304 bp); mean depth
      0.012 and 220. Four near-empty pools (KS-2-83, KS-2-90, SN-1-82,
      SN-2-91) have outputs under 1 MB. Copied to
      `/Volumes/lan_images4/tany_time/bam-cov/` (62 files, `gzip -t` OK);
      `_data/bed-summs.rds` rebuilt from them (`_R/sample-qc.R`): median
      depth 0-231 (time series coverage-30 pools 36-163). In the pools
      with the most mitochondrial reads, ~150-270 bp at
      CONTIG37:2,739,348 reach 0.66-2.2 million reads (SN-2-96,
      KS-1-15-B); `samtools mpileup` caps depth at 8,000 and `01-snape`
      masks each pool's sites above its 95th coverage percentile.
- [x] **G4.** *(parallel with G2)* `05-kraken` (1 job, all pools) →
      `./check-step.sh 01-variant-calling/05-kraken` (pools classified so
      far, warnings). Run 2026-10-30 (cluster 12391035): PASS, 1.2 h, all
      62 pools classified with both databases. Summarized with
      `_R/kraken-summary.R` → `_data/kraken-*.csv`,
      `_plots/kraken-groups.pdf`. The SILVA host denominator is empty (see
      `TODO.md`, Parasites with Kraken2); replaced by reads per million
      unmapped and per million total reads.

## H. Allele frequencies (`02-allele-freqs/`)

- [ ] **H1.** `01-snape` (62 jobs) → `./check-step.sh 02-allele-freqs/01-snape`
      (also lists contigs where SNAPE reported errors).
- [ ] **H2.** `02-npstat` (62 jobs; after H1, since it reads each
      pool's `snape/<pool>_snape.tar`; can run alongside H3–H6) →
      `./check-step.sh 02-allele-freqs/02-npstat`.
- [ ] **H3.** `03-pool-qc` (1 job) → `./check-step.sh 02-allele-freqs/03-pool-qc`.
      Since 2026-10-30 it reads every candidate pool's whole BAM to count
      reads on the CONTIG37 stretch, so it runs longer than before.
- [ ] **H4.** Locally: download `$S/dna/pool-qc/pool-qc.tar` and extract it
      to `pool-qc/` in the local data folder (`_R/data-dir.R`; see
      `CLAUDE.md`); run `Rscript _R/pool-qc.R` (writes
      `_bash/samples/{time,space}-pools.txt` by the pre-set outlier rule)
      and `Rscript _R/mito-contam.R`. After G4, extract
      `$S/dna/kraken/kraken.tar` to `kraken/` in the local data folder and
      run `Rscript _R/kraken-summary.R`. Commit and push the pool lists.
- [ ] **H5.** On the submit server, `git pull` (gets the pool lists).
- [ ] **H6.** `04-combine-syncs` (1 job) →
      `./check-step.sh 02-allele-freqs/04-combine-syncs`.

## I. BayPass (`03-select-time/`)

- [ ] **I1.** `00-baypass-inputs` (3 jobs) →
      `./check-step.sh 03-select-time/00-baypass-inputs` (fails if a
      setting keeps no SNPs).
- [ ] **I2.** `01-baypass-init` (6 jobs, 25–50 h each; one first) →
      `./check-step.sh 03-select-time/01-baypass-init`.
- [ ] **I3.** `02-check-omegas` (1 job) →
      `./check-step.sh 03-select-time/02-check-omegas` (PASS or WARN on Ω
      consistency; covariate runs use the all-SNP Ω either way).
- [ ] **I4.** `03-baypass-beta` (60 jobs, 25–50 h each; one first) →
      `./check-step.sh 03-select-time/03-baypass-beta`.

## J. After the rerun

- [ ] **J1.** Update `_R/baypass-analysis.R`, `baypass-check-omegas.R`,
      and `baypass-check-candidates.R` for the new output names and the
      missing-data grid.
- [ ] **J2.** Continue with the review steps in `TODO.md`.
