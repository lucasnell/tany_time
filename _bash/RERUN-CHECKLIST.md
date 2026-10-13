# Cluster rerun checklist

Full rerun of the pipeline with container v2.2.0. Paths are on CHTC
staging (`/staging/l/lnell/`, written `$S` below) unless noted. Submit each
step from inside its folder after `mkdir -p logs/<step>`. Details for each
step are in that folder's `README.md`.

Steps marked *(parallel)* can run alongside the step before them.

## A. Local preparation

- [ ] **A1.** Pull the latest repo (`git pull`); check `git status` is clean.
- [ ] **A2.** Regenerate job lists (only needed if samples, QC, or the
      BayPass grid changed): `Rscript _R/make-job-lists.R`.
- [ ] **A3.** Make the sample-info tar for BayPass and copy it to
      staging (commands at the top of
      `03-select-time/00-baypass-inputs.sh`) → `$S/tany-sample-info.tar.gz`.

## B. Staging setup (on the transfer/submit server)

- [ ] **B1.** Create output folders (only pool QC, Kraken2, and the SILVA
      build create their own):

      ```bash
      S=/staging/l/lnell
      mkdir -p $S/dna/{raw_fq,trimmed,bwa,mpileup,ont,snape,npstat,pool-qc,kraken} \
          $S/dna/baypass/{init,beta} $S/{other-genomes,silva,kraken2-db}
      ```

- [ ] **B2.** Confirm inputs already on staging:
  - [ ] Raw reads: `$S/dna/raw_fq/<pool>_L002_R{1,2}_001.fastq.gz` for all
        62 pools in `samples/all-pools.txt`.
  - [ ] Nanopore reads: `$S/dna/ont/basecalls_guppy-5.0.11.fastq.gz`.
  - [ ] Assembly: `$S/Tgraci_assembly.fasta.gz`.
  - [ ] Repeats: `$S/Tgraci_repeats_locs.gff3.gz`.
- [ ] **B3.** Copy the other assemblies from
      `~/_data/midgenomes/assemblies/` to `$S/other-genomes/`:
      `Cripar_assembly.fasta.gz`, `Ctenta2023_assembly.fasta.gz`,
      `Saterr_assembly.fasta.gz`, `Pakamu_assembly.fasta.gz`.
- [ ] **B4.** Download the three SILVA 138.2 files into `$S/silva/` (URLs at
      the top of `01-variant-calling/06a-silva-db.sh`).
- [ ] **B5.** Download Kraken2 PlusPF-16 into `$S/kraken2-db/` from
      https://benlangmead.github.io/aws-indexes/k2, and set
      `KRAKEN_DB_TAR` in `01-variant-calling/06b-kraken.sh` to its file
      name (currently a placeholder, `k2_pluspf_16gb.tar.gz`).
- [ ] **B6.** Copy the whole `_bash/` folder to the submit server, keeping
      its layout (`.sub` files read `../samples/` and `../references/`).

## C. Container

- [ ] **C1.** Build the image (`_apptainer/README.md`): copy
      `_apptainer/` files to the submit server, run `_build.sh` in an
      interactive build job → `$S/tany_time.sif`.
- [ ] **C2.** Quick check in an interactive job: each environment
      activates (`main-env`, `qc-env`, `gatk-env`, `gatk3-env`,
      `baypass-env`, `kraken-env`, `mito-env`), `/app/bin` scripts are on
      `PATH`, and `medaka tools list_models` (in `mito-env`) shows
      `r941_min_hac_g507`.

## D. References (`01-variant-calling/`)

- [ ] **D1.** `00a-mitogenome` (1 job).
- [ ] **D2.** **Check** `coi-check.tsv` in `$S/Tgraci_mitogenome.tar`:
      both T. gracilentus COI barcodes align over their full length at
      ~100% identity, and the mitogenome is ~16 kb. Stop here if not.
- [ ] **D3.** `00b-hologenome` (1 job) → `$S/hologenome.tar`. Its log
      lists sequences and bp per group; check every species is there.
- [ ] **D4.** *(parallel with D1–D3)* `06a-silva-db` (1 job) →
      `$S/kraken2-db/silva-138.2.tar.gz`.
- [ ] **D5.** *(optional, parallel)* `04-ont-pi` (1 job). The SNAPE prior
      θ = 0.0088 is hard-coded in `02-allele-freqs/01-snape.sh`; only
      update it if this gives a clearly different value.

## E. Reads to pileups (`01-variant-calling/`, 62 jobs each)

- [ ] **E1.** `01-fastp` *(can start during D)*.
- [ ] **E2.** `02-bwa` (needs `hologenome.tar` from D3). Spot-check a few
      `<pool>_species-counts.tsv` in `$S/dna/bwa/<pool>_bwa.tar.gz`.
- [ ] **E3.** `03-mpileup` (9+ h per job). Check the first job finishes
      well within the 72-h limit.
- [ ] **E4.** *(parallel with E3)* `05-bam-cov`.
- [ ] **E5.** *(parallel with E3)* `06b-kraken` (1 job, all pools; long:
      check its log after a few pools).

## F. Allele frequencies (`02-allele-freqs/`)

- [ ] **F1.** `01-snape` (62 jobs).
- [ ] **F2.** *(parallel)* `02-npstat` (62 jobs).
- [ ] **F3.** `03-pool-qc` (1 job) → `$S/dna/pool-qc/pool-qc.tar`.
- [ ] **F4.** Locally: download `pool-qc.tar` and extract it to
      `~/_data/_tany-time/pool-qc/`; run `Rscript _R/pool-qc.R` (needs
      `rrcov`). It writes `_bash/samples/{time,space}-pools.txt`. Review
      `_plots/pool-qc-*.pdf` and the excluded pools, then commit the lists.
- [ ] **F5.** Locally: `Rscript _R/mito-contam.R` (contamination index;
      uses the same `pool-qc` folder). After E5, download `kraken.tar`
      to `~/_data/_tany-time/kraken/` and run `Rscript _R/kraken-summary.R`;
      read its warning about SILVA names that matched nothing.
- [ ] **F6.** Copy the new `samples/{time,space}-pools.txt` to the submit
      server.
- [ ] **F7.** `04-combine-syncs` (1 job).

## G. BayPass (`03-select-time/`)

- [ ] **G1.** `00-baypass-inputs` (3 jobs, one per missing-data setting).
      Check each log's "SNPs read / kept" counts.
- [ ] **G2.** `01-baypass-init` (6 jobs). Submit one first and check its
      progress against the 72-h limit (expected 25–50 h).
- [ ] **G3.** Check Ω consistency for `mp075` (subsamples 1–3 vs. all
      SNPs; FMD distance) before going on. `_R/baypass-check-omegas.R`
      needs updating for the new file names first.
- [ ] **G4.** `02-baypass-beta` (60 jobs). Submit one first and check its
      progress, then the rest.

## H. After the rerun

- [ ] **H1.** Update `_R/baypass-analysis.R`, `baypass-check-omegas.R`,
      and `baypass-check-candidates.R` for the new output names and the
      missing-data grid.
- [ ] **H2.** Continue with the review steps in `TODO.md`.
