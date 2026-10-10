# `_bash/`: cluster pipeline

Scripts for the UW–Madison CHTC cluster (HTCondor). Each step is a `.sh`
script plus a `.sub` submit file, numbered in run order within each folder.
All jobs run inside the Apptainer container built from `_apptainer/`.

| Folder | Contents |
|---|---|
| `_apptainer/` | Container definition, conda environments, shared shell functions, and shared scripts |
| `samples/` | Job lists (pool names, pool sizes, QC-passing pools) |
| `01-variant-calling/` | Reads → trimmed reads → BAM → mpileup |
| `02-allele-freqs/` | mpileup → SNAPE allele frequencies, counts, and diversity statistics |
| `03-select-time/` | BayPass: allele frequencies vs. abundance through time |

## Run order

1. Build the container: `_apptainer/_build.sh` → `/staging/lnell/tany_time.sif`.
2. Write job lists locally: `Rscript _R/make-job-lists.R` (from the repo root).
3. `01-variant-calling/`: `01-fastp` → `02-bwa` → `03-mpileup`;
   `05-bam-cov` after `02-bwa`. `04-ont-pi` is independent (one job on
   Nanopore reads) and only sets the SNAPE prior θ used in `02-allele-freqs/01-snape`.
4. `02-allele-freqs/`: `01-snape` → `03-combine-syncs`; `02-npstat` after `01-snape`.
5. `03-select-time/`: `00-baypass-inputs` → `01-baypass-init` → `02-baypass-beta`.

Each folder's `README.md` lists inputs, outputs, and parameters.

## Running on CHTC

- Copy the whole `_bash/` folder to the submit server, keeping its layout:
  `.sub` files read job lists from `../samples/`.
- Submit from inside a step's folder, after creating its log folder, e.g.:

  ```bash
  cd 01-variant-calling
  mkdir -p logs/01-fastp
  condor_submit 01-fastp.sub
  ```

- Inputs and outputs are on `/staging/lnell/` (paths are set near the top of
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
