# `_apptainer/`: container

`tany_time.def` builds `tany_time.sif` (v2.3.0; v2.2.0 had medaka model `r941_min_hac_g507` instead) from
`condaforge/miniforge3:26.7.2-0` (Ubuntu 24.04).

## Build

```bash
# From the cluster (after updating via `git pull`):
cd /home/lnell/tany_time/_bash/_apptainer
mkdir -p /staging/l/lnell/_apptainer/
cp -r *.yml helpers.sh *.def bin _build.sh /staging/l/lnell/_apptainer/
# On the cluster, in an interactive build job:
cd /staging/l/lnell/_apptainer && chmod +x ./_build.sh && ./_build.sh     # moves tany_time.sif to /staging/l/lnell/
```

## Conda environments

Versions are pinned in each `.yml` (dry-run solved for linux-64 with strict
channel priority; channels conda-forge then bioconda).

| Environment | Main tools | Used by |
|---|---|---|
| `main-env` | fastp, bbmap, bwa, samtools, bcftools, bedtools, minimap2, vcftools, python (numpy, pandas) | most steps |
| `qc-env` | fastqc, kat, mapdamage2 | `01-fastp` |
| `gatk-env` | GATK 4 (Picard tools) | `02-bwa`, `03-mpileup` |
| `gatk3-env` | GATK 3.8 (IndelRealigner) | `03-mpileup` |
| `baypass-env` | BayPass 3.1, R with poolfstat and tidyverse | `03-select-time/` |
| `kraken-env` | Kraken2 2.17.2, samtools, seqtk | `00-references/03-silva-db`, `01-variant-calling/05-kraken` |
| `mito-env` | Flye 2.9.6, medaka 2.2.2 (R9.4.1 model `r941_min_sup_g507` downloaded at build), minimap2, bwa, polypolish 0.7.1, samtools, seqtk | `00-references/01-mitogenome` |

Also built from source: `snape-pooled` (DEST-bio fork, commit `284bfe0`;
the original repo was deleted) and `npstat` (commit `dc230fb`).

## Shared shell functions (`helpers.sh`)

Appended to `/app/.bashrc`, so scripts get them via `. /app/.bashrc`.

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

## Shared scripts (`bin/`, on `PATH` as `/app/bin`)

| Script | Purpose | Used by |
|---|---|---|
| `PickleRef.py`, `Mpileup2Sync.py`, `SNAPE2SYNC.py`, `MaskSYNC_snape.py` | SNAPE helper scripts from DEST (Kapun et al. 2021); see `02-allele-freqs/README.md` for local edits | `01-snape` |
| `combine-syncs.py` | Combine per-pool gSYNC files; counts by default, SNAPE info with `--snape`; `--min-pp` keeps sites polymorphic in at least one pool; `--clean-mono` assigns reads to the reference where SNAPE calls a pool monomorphic | `04-combine-syncs` |
| `pool-qc-sync.py` | Per-pool missing-data rate, SNPs, and private SNPs from partially masked gSYNC files | `03-pool-qc` |
| `pool-qc-files.py` | Per-pool GC, duplication, and mapping rates (including reads assigned to *T. gracilentus* vs. other species) from fastp, Picard, bamtools, and `02-bwa` outputs | `03-pool-qc` |
| `mito-collapse.sh` | Collapse the assembled mitochondrial contig to the gene region, the unique control-region sequence on each side, and one repeat unit (needs minimap2, samtools: `mito-env`). The job transfers this file from the repo and runs it, so the container's copy (v2.3.0, before the 2026-10-27 fix) is unused | `00-references/01-mitogenome` |
| `run-parallel.py` | Run a command once per item (e.g., contig) in parallel and report failures | `01-snape`, `02-npstat` |

Python scripts run with whichever environment is active (`main-env` in all
current uses).
