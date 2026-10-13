# `_bash/`: cluster pipeline

Scripts for the UW–Madison CHTC cluster (HTCondor). Each step is a `.sh`
script plus a `.sub` submit file, numbered in run order within each folder.
All jobs run inside the Apptainer container built from `_apptainer/`.

| Folder | Contents |
|---|---|
| `_apptainer/` | Container definition, conda environments, shared shell functions, and shared scripts |
| `samples/` | Job lists (pool names, pool sizes, QC-passing pools) |
| `references/` | Small reference sequences (mitochondrial bait, COI barcodes) |
| `tests/` | Local checks of pipeline behavior (results in `_data/`) |
| `00-references/` | References built once: *T. gracilentus* mitogenome, combined reference for competitive mapping, SILVA Kraken2 database, nanopore θ estimate |
| `01-variant-calling/` | Reads → trimmed reads → BAM → mpileup |
| `02-allele-freqs/` | mpileup → SNAPE allele frequencies, counts, and diversity statistics |
| `03-select-time/` | BayPass: allele frequencies vs. abundance through time |

## Run order

Step-by-step checklist for a full rerun: `RERUN-CHECKLIST.md`.

1. Build the container: `_apptainer/_build.sh` → `/staging/l/lnell/tany_time.sif`.
2. Write job lists locally: `Rscript _R/make-job-lists.R` (from the repo root).
3. `01-variant-calling/01-fastp` (the mitogenome is polished with trimmed
   reads from one pool).
4. `00-references/`: `01-mitogenome` → `02-hologenome`; `03-silva-db` and
   `04-ont-pi` (optional) are independent.
5. `01-variant-calling/`: `02-bwa` → `03-mpileup`; `04-bam-cov` and
   `05-kraken` after `02-bwa`.
6. `02-allele-freqs/`: `01-snape` → `03-pool-qc` → locally, `_R/pool-qc.R`
   (choose pools from QC metrics; writes `samples/{time,space}-pools.txt`)
   → `04-combine-syncs`. `02-npstat` runs any time after `01-snape`.
7. `03-select-time/`: `00-baypass-inputs` → `01-baypass-init` →
   `02-check-omegas` → `03-baypass-beta`.

After (or while) a step runs, check it with `check-step.sh` from the
submit server, e.g. `./check-step.sh 01-variant-calling/02-bwa`. It
reports finished, failed, running, and unsubmitted jobs; the longest run
time (warning above 54 h; CHTC's limit is 72 h); outputs missing from
staging; and a step-specific summary. It ends with `RESULT: PASS` or the
number of problems.

Each folder's `README.md` lists inputs, outputs, and parameters.

## Running on CHTC

- The repo is cloned on the submit server at `/home/lnell/tany_time`
  (public, so HTTPS works). Jobs are submitted from inside the clone, so
  `.sub` files find `../samples/` and `../references/`:

  ```bash
  # once:
  cd /home/lnell && git clone https://github.com/lucasnell/tany_time.git
  # before each session (after pushing local changes):
  cd /home/lnell/tany_time && git pull
  ```

- Submit from inside a step's folder, after creating its log folder, e.g.:

  ```bash
  cd /home/lnell/tany_time/_bash/01-variant-calling
  mkdir -p logs/01-fastp
  condor_submit 01-fastp.sub
  ```

- `logs/` folders are git-ignored. Edit files locally and push, not on
  the cluster, so `git pull` never conflicts. Job scripts are stored as
  executable in git.

- Inputs and outputs are on `/staging/l/lnell/` (paths are set near the top of
  each `.sh` under "Edit these for your system"). Steps that set
  `output_destination` in their `.sub` also copy outputs to ResearchDrive
  (`pelican://chtc.wisc.edu/researchdrive/arives/CHTC/`).
- CHTC's default run limit is 72 h per job. The longest steps are
  `03-mpileup` (indel realignment, 9+ h) and BayPass (~25–50 h with
  `NVAL=5000`; see `03-select-time/README.md`).
- Every `.sh` starts with `. /app/.bashrc` (loads the functions in
  `_apptainer/helpers.sh`) and `conda activate <env>`.

## Conventions

- `READ_BASE` / `readbase` is the pool name used by the sequencing center
  (e.g., `SN-2-07_S31`); raw reads are
  `${READ_BASE}_L002_R{1,2}_001.fastq.gz`.
- Scripts check inputs with `require_files` and check each command with
  `check_exit_status`, which prints logs and exits on failure.
- Scripts that copy outputs to ResearchDrive end with `exit 0` rather than
  `safe_exit 0`, because `safe_exit` deletes the scratch directory before
  HTCondor can transfer outputs.
