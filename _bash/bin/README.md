# `bin/`: shared scripts and shell functions

Every job transfers this folder from the repo (`transfer_input_files =
..., ../bin` in each `.sub`), so edits take effect after `git pull` on the
submit server, with no container rebuild. The container (`_apptainer/`)
holds only software. Container images up to v2.4.0 also had copies in
`/app/bin` and `/app/.bashrc`; jobs ignore those because `helpers.sh` is
sourced after `/app/.bashrc` and puts this folder first on `PATH`.

## Shell functions (`helpers.sh`)

Each job sources it after `. /app/.bashrc` (`. bin/helpers.sh`). Sourcing
it also puts this folder first on `PATH`.

| Function | Purpose |
|---|---|
| `count_threads` | Number of CPUs assigned by HTCondor |
| `safe_exit` | Exit; in batch jobs also clears the scratch directory |
| `check_exit_status` | Exit with logs if the previous command failed |
| `require_files` | Exit if any listed input file is missing |
| `stage_genome` | Copy and decompress a FASTA (`upper` removes softmasking) |
| `contig_names` | Sequence names from a FASTA, in order |
| `split_by_contig` | Split a file into one file per contig, checking names against `CONTIG_NAMES` |
| `read_tar_name` | Name of read file 1 or 2 inside a tar |
| `call_bam_stats` | `bamtools stats` on a BAM (needs `main-env`) |

## Scripts

| Script | Purpose | Used by |
|---|---|---|
| `PickleRef.py`, `Mpileup2Sync.py`, `SNAPE2SYNC.py`, `MaskSYNC_snape.py` | SNAPE helper scripts from DEST (Kapun et al. 2021); see `02-allele-freqs/README.md` for local edits | `01-snape` |
| `combine-syncs.py` | Combine per-pool gSYNC files; counts by default, SNAPE info with `--snape`; `--min-pp` keeps sites polymorphic in at least one pool; `--clean-mono` assigns reads to the reference where SNAPE calls a pool monomorphic | `04-combine-syncs` |
| `pool-qc-sync.py` | Per-pool missing-data rate, SNPs, and private SNPs from partially masked gSYNC files | `03-pool-qc` |
| `pool-qc-files.py` | Per-pool GC, duplication, and mapping rates (including reads assigned to *T. gracilentus* vs. other species) from fastp, Picard, bamtools, and `02-bwa` outputs | `03-pool-qc` |
| `mito-collapse.sh` | Collapse the assembled mitochondrial contig to the gene region, the unique control-region sequence on each side, and one repeat unit (needs minimap2, samtools: `mito-env`). | `00-references/01-mitogenome` |
| `run-parallel.py` | Run a command once per item (e.g., contig) in parallel and report failures | `01-snape`, `02-npstat` |

Python scripts run with whichever environment is active (`main-env` in all
current uses).
