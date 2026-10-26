# Cluster rerun checklist

Full rerun of the pipeline with container v2.3.0 (v2.2.0 plus the sup
medaka model; built v2.2.0 for steps run so far). Paths are on CHTC
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

- [ ] **C3.** Rebuild as v2.3.0 (medaka model `r941_min_sup_g507` for the
      sup reads, 2026-10-23; `bin/mito-collapse.sh`, 2026-10-26) as in C1,
      then C2 with `grep r941_min_sup_g507` and `which mito-collapse.sh`.
      Needed before F1.

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
For steps marked "one first", submit a single job (e.g., by temporarily
keeping only one line of the job list) and check its run time before
submitting the rest.

## E. Trimming (`01-variant-calling/`)

- [x] **E1.** `01-fastp` (62 jobs) → `./check-step.sh 01-variant-calling/01-fastp`.

## F. References (`00-references/`)

- [ ] **F1.** `01-mitogenome` (1 job; needs trimmed reads of
      SN-1-08_S24 from E1). It fails by itself, writing nothing, unless
      COI barcode MZ659759.1 matches the collapsed mitogenome
      (gene region, unique control-region sides, one repeat unit) →
      `./check-step.sh 00-references/01-mitogenome`.
- [ ] **F2.** `02-hologenome` (1 job). It fails by itself unless every
      species, the mitogenome, and every barcode group are present →
      `./check-step.sh 00-references/02-hologenome`.
- [x] **F3.** *(parallel with F1–F2)* `03-silva-db` (1 job) →
      `./check-step.sh 00-references/03-silva-db`.
- [x] **F4.** *(optional, parallel)* `04-ont-pi` (1 job) →
      `./check-step.sh 00-references/04-ont-pi` compares the new θ with
      the SNAPE prior in `02-allele-freqs/01-snape.sh` and warns if they
      differ by more than 10%.

## G. Reads to pileups (`01-variant-calling/`)

- [ ] **G1.** `02-bwa` (62 jobs) → `./check-step.sh 01-variant-calling/02-bwa`
      (also tabulates each pool's reads on *T. gracilentus*, other
      species, the mitogenome, and barcodes).
- [ ] **G2.** `03-mpileup` (62 jobs, 9+ h each; one first) →
      `./check-step.sh 01-variant-calling/03-mpileup`.
- [ ] **G3.** *(parallel with G2)* `04-bam-cov` (62 jobs) →
      `./check-step.sh 01-variant-calling/04-bam-cov`.
- [ ] **G4.** *(parallel with G2)* `05-kraken` (1 job, all pools) →
      `./check-step.sh 01-variant-calling/05-kraken` (pools classified so
      far, warnings).

## H. Allele frequencies (`02-allele-freqs/`)

- [ ] **H1.** `01-snape` (62 jobs) → `./check-step.sh 02-allele-freqs/01-snape`
      (also lists contigs where SNAPE reported errors).
- [ ] **H2.** *(parallel)* `02-npstat` (62 jobs) →
      `./check-step.sh 02-allele-freqs/02-npstat`.
- [ ] **H3.** `03-pool-qc` (1 job) → `./check-step.sh 02-allele-freqs/03-pool-qc`.
- [ ] **H4.** Locally: download `$S/dna/pool-qc/pool-qc.tar` and extract it
      to `~/_data/_tany-time/pool-qc/`; run `Rscript _R/pool-qc.R` (writes
      `_bash/samples/{time,space}-pools.txt` by the pre-set outlier rule)
      and `Rscript _R/mito-contam.R`. After G4, download
      `$S/dna/kraken/kraken.tar` to `~/_data/_tany-time/kraken/` and run
      `Rscript _R/kraken-summary.R`. Commit and push the pool lists.
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
