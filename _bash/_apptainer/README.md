# `_apptainer/`: container

`tany_time.def` builds `tany_time.sif` (v2.5.0, not yet built: project
scripts and shell functions removed from the image; v2.4.0: bbmap back to 39.28, since
40.02's BBMerge loses or hangs on paired reads; v2.3.0 had bbmap 40.02;
v2.2.0 had medaka model `r941_min_hac_g507`) from
`condaforge/miniforge3:26.7.2-0` (Ubuntu 24.04).

## Build

```bash
# From the cluster (after updating via `git pull`):
cd /home/lnell/tany_time/_bash/_apptainer
mkdir -p /staging/l/lnell/_apptainer/
cp -r *.yml *.def _build.sh /staging/l/lnell/_apptainer/
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

Project scripts and shell functions are not in the image; they are in
`_bash/bin/` (see its `README.md`), and each job transfers them from the
repo, so editing them needs no rebuild.
