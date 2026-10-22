# `00-basecall/`: re-basecall the nanopore reads

Re-basecalls the raw nanopore signal of the genome individual (R9.4.1
flow cell, sequenced 2019; first basecalled with Guppy 5.0.11 hac) with
Dorado 0.9.6 and the R9.4.1 sup model. Dorado 0.9.6 is the last release
that reads FAST5 and has R9.4.1 models (both removed in Dorado 1.0.0). On a
local test of 200 reads, sup's mean predicted quality was 17.3 vs 14.4 for
Guppy hac (`../tests/dorado-r941-test.sh` → `_data/dorado-r941-test.txt`).
`../check-step.sh 00-basecall/<step>` summarizes any step.

| Step | Does | Input (`/staging/l/lnell/`) | Output (`/staging/l/lnell/`) |
|---|---|---|---|
| `00-setup` | Run once in an interactive build job (`condor_submit -i 00-setup.sub`, then `bash 00-setup.sh`): download Dorado 0.9.6 (linux-x64) and model `dna_r9.4.1_e8_sup@v3.6` | (internet) | `dorado/` |
| `01-split-fast5` | Stream the archive and write uncompressed tars of 20 FAST5 files (~8 GB) each; checks the archive's MD5 in the same pass. One CPU job, several hours | `dna/ont/ives_fast5.tar.bz2` | `dna/ont/fast5_chunks/chunk_NNN.tar`, `manifest.tsv` |
| `02-dorado` | Basecall one chunk per GPU job (CHTC GPU Lab); runs without the container | `dorado/`, `dna/ont/fast5_chunks/` | `dna/ont/dorado/chunk_NNN.bam`, `chunk_NNN_summary.tsv` |
| `03-combine` | Reads with mean quality >= 10 to one FASTQ; all reads to one summary table | `dna/ont/dorado/` | `dna/ont/basecalls_dorado-0.9.6_sup.fastq.gz`, `basecalls_dorado-0.9.6_sup--summary.tsv.gz` |

Before `01-split-fast5`:

1. Check staging space (archive 330 GB, chunks about the same, BAMs tens
   of GB); the archive can be deleted once `01-split-fast5` passes its MD5
   check.
2. Copy the archive from the external drive
   (`/Volumes/TanyGen5TB/nanopore/ives_fast5.tar.bz2`) to
   `/staging/l/lnell/dna/ont/`, e.g. with `rsync --partial --progress`
   through CHTC's transfer server.
3. Start an interactive build job (`mkdir -p logs/00-setup && condor_submit -i
   00-setup.sub`) and run `bash 00-setup.sh` inside it, then `exit`.

Submit `02-dorado` after `01-split-fast5` finishes: its submit file queues
one job per `chunk_*.tar` on staging. The GPU settings in `02-dorado.sub`
(`+WantGPULab`, `+GPUJobLength`, `gpus_minimum_capability`,
`gpus_minimum_memory`) follow CHTC's GPU guide
(https://chtc.cs.wisc.edu/uw-research-computing/gpu-jobs, checked
2026-10-22).

Not yet done: the pipeline steps that use nanopore reads (`01-mitogenome`,
`04-ont-pi` in `00-references/`) still read the Guppy basecalls. Switching
them to the sup basecalls also means switching medaka to a sup model.
