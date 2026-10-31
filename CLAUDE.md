# tany_time

Pool-seq population genomics of the midge *Tanytarsus gracilentus* at Mývatn,
Iceland. Two sample sets:

- **Time series**: archived emergence-trap samples from site Syðri Neslönd
  (SN), 17 pools spanning 1992–2015, two generations per year (`a` = before
  mid-July, `b` = after).
- **Space**: 18 pools from Mývatn sites and other Icelandic lakes, 2019.

Goal: understand how population crashes affect *T. gracilentus* adaptive
evolution and neutral genetic variation. Crashes specifically are the main
interest; association with abundance in general is secondary.

Current analysis: SNPs whose allele frequencies covary with abundance
(log(1 + N)) across the time series, tested with BayPass, then matched to
genes, BLASTed to *D. melanogaster*, and passed to a time-series analysis
aimed at balancing selection. A 2026-10 review found this candidate set is
not yet interpretable (see `TODO.md`).

## Study facts (from the owner)

- Abundance (`tany-abundances.csv`) is from the SN trap alone.
- Samples are males only; species identity was checked by male genitalia.
  No other *Tanytarsus* species live in Mývatn.
- Deep-crash generations (e.g. 1984, 2003–2005, 2011) had too few adults to
  sequence; no further pools can be added, and existing pools can't be
  re-sequenced (no material left).
- Midge community: 37 chironomid species at Mývatn; trap counts in
  Einarsson et al. 2026, *Aquatic Ecology* 60:58 (Table 1;
  doi:10.1007/s10452-026-10294-2; Zotero `WEMS8F4U`). Most abundant in
  traps: *T. gracilentus*, *Cricotopus sylvestris*, *Simulium vittatum*
  (black fly), *Eukiefferiella minor*, *Micropsectra lindrothi*,
  *Procladius islandicus*, *Chironomus islandicus*.
- Archived samples were stored in plastic containers in a basement for
  years; some were stored poorly (some dried out), which likely degraded
  DNA. Older samples were probably more affected, but quality doesn't track
  year reliably.

See `PROJECT_INDEX.md` for the file map and `TODO.md` for open items.

## Layout

- `_bash/` — cluster pipeline (UW–Madison CHTC, HTCondor); see
  `_bash/README.md`. Each step is a `.sh` + `.sub` pair, numbered in run
  order. Runs inside the Apptainer image built from
  `_bash/_apptainer/tany_time.def` (`/staging/l/lnell/tany_time.sif`).
  Shared shell functions (`helpers.sh`) and scripts are in `_bash/bin/`,
  not in the container: every job transfers that folder and sources
  `bin/helpers.sh`, which puts it first on `PATH` (owner's rule,
  2026-10-31: one copy of each script, in the repo). Small reference sequences for jobs are in `_bash/references/`
  (from `_R/get-mito-refs.R`). Job lists in `_bash/samples/` come from `_R/make-job-lists.R`,
  except the included-pool lists (`time-pools.txt`, `space-pools.txt`),
  which come from `_R/pool-qc.R`.
  `_bash/tests/` holds local checks of pipeline behavior.
  `_bash/00-basecall/` re-basecalls the raw nanopore signal (archive
  `/Volumes/TanyGen5TB/nanopore/ives_fast5.tar.bz2`) with Dorado 0.9.6, the
  last release with R9.4.1 models and FAST5 input; `02-dorado` runs on GPUs
  without the container.
  CHTC disables HTCondor email notification (submissions with
  `notification`/`notify_user` are rejected; 2026-10-23).
  `_bash/check-step.sh <folder>/<step>` (on the submit server) checks a
  step's jobs, run times, outputs, and step-specific results; checks
  should be scripted like this rather than left to manual inspection
  (owner's request, 2026-10-13).
  Software building or installing (container builds, downloads of tools or
  models) runs in an interactive build job (`condor_submit -i` with
  `+IsBuildJob = true`; https://chtc.cs.wisc.edu/uw-research-computing/inter-submit),
  never on the submit server (owner's rule, 2026-10-19).
  `_bash/zzz-old/` (git-ignored, local only) holds old assembly scripts.
  Conda envs: `main-env`, `qc-env`, `gatk-env`, `gatk3-env`, `baypass-env`,
  `kraken-env`, `mito-env`.
  Pins updated to v2.0.0 on 2026-10-09 (BayPass 3.1, fastp 1.4.0, samtools
  1.24, …); all existing outputs were made with v1.0.0 (BayPass 2.31).
  bbmap stays at 39.28: 40.02's BBMerge silently drops or hangs on paired
  reads (container v2.4.0, 2026-10-29). Check tool outputs, not just exit
  status, when changing pins.
  snape-pooled is built from the DEST-bio fork (original repo deleted).
- `_R/` — current (2025) local R analysis. Run from the repo root.
- `_data/` — small inputs and derived datasets. `*.rds`, `*.gz`, and `*.pdf`
  are git-ignored, so derived `.rds` files exist only locally.
- `_plots/` — figures (git-ignored PDFs).
- `_refs/` — reference documents (BayPass v2.31 and v3.1 manuals; git-ignored
  PDFs).
- `_legacy/` — older analyses (mostly 2020–2022) formerly at the repo root (previous
  MAPGD/SNAPE pipeline, spatial FST, admixture, simulations). They read from
  old paths (`~/Box Sync`, `~/Stanford_Drive`, `~/Cornell_Drive`,
  `full-DNA-info.csv`), and their relative paths (`sourceCpp("helpers.cpp")`,
  `_plots/`, `.Rprofile`) assumed the repo root, so they will not run as-is.
  Do not edit unless asked. `_R/crash-test.R` is also legacy (superseded by
  BayPass).

## Data outside the repo

Scripts read large files by absolute path; none are in git.

**Local data folder** (owner's rule, 2026-10-30): cluster outputs that
local analyses need are copied, as soon as they exist, into
`/Volumes/lan_images4/tany_time/<output>/` (e.g., `bam-cov/`, `kraken/`,
`pool-qc/`); if only `LAN_4TB_SSD` is connected, into
`/Volumes/LAN_4TB_SSD/tany_time/` instead. Whenever both drives are
connected, sync the two folders both ways so each holds everything (both
exFAT; add, never delete):
`rsync -rt --update --modify-window=1 /Volumes/lan_images4/tany_time/ /Volumes/LAN_4TB_SSD/tany_time/`,
then the same with the two paths swapped. R scripts find the folder with
`_R/data-dir.R` (`tany_dir()`). Outputs from before this rule (SNAPE,
BLAST, trimmed reports) stay in `~/_data/_tany-time/` until rerun.

| Path | Contents |
|---|---|
| `~/_data/_tany-time/snape/{snape,count}_masked_{time,space,all}/` | Combined SNAPE frequency and read-count sync files |
| `~/_data/_baypass/{baypass-inputs,init,beta}/` | BayPass inputs and outputs (local copy holds only `beta0` runs 1–10) |
| `~/_data/midgenomes/` | Assemblies (`assemblies/<code>_assembly.fasta.gz`; `Ctenta2023` and `Saterr` added 2026-10-12 from NCBI, GCA_963573255.1 and GCA_033063855.1), GO terms, OrthoFinder HOG→GO tables |
| `~/_data/Tgraci_genes.gff3` | Gene annotation |
| `~/_data/_tany-time/blast/` | Protein FASTA (`Tgraci_proteins.faa`), candidate-gene proteins (`Tgraci-funct.faa`), and BLASTP output (`Tgraci-funct-Dmel.csv`) used in `baypass-analysis.R` and `make-time-series-plot.R` (moved from `~/Desktop` 2026-10-13); BLAST database in `~/blastdb/` |
| `~/_data/_tany-time/trimmed/trimmed_reports/` | fastp/FastQC per-sample reports and `multiqc_data/` (`sample-qc.R` still points to the old `~/_data/trimmed_reports`) |
| Local data folder (`tany_time/` on lan_images4 or LAN_4TB_SSD; above) | `bam-cov/`: per-pool coverage on *T. gracilentus* contigs (`<pool>_bwa_cov.bed.gz`, from `04-bam-cov`, 2026-10-30 run; read by `_R/sample-qc.R`; older coverage from the previous mapping is on `/Volumes/Lucas512`). `kraken/`: Kraken2 reports from `05-kraken` (`_R/kraken-summary.R`). `pool-qc/`: `03-pool-qc` output (`_R/pool-qc.R`, `_R/mito-contam.R`) |
| `/Volumes/lan_images4/` (external drive, exFAT) | `ont/`: Dorado sup basecalls and summary (`00-basecall/03-combine`); `assembly/{flye,necat,nextdenovo,smartdenovo}/contigs_*.fasta.gz`: 2019 whole-genome assemblies (copied from `/Volumes/TanyGen5TB/assembly/` 2026-10-23, MD5s match); `ill/dna/trimmed/` and `rna/trimmed/`: fastp-trimmed DNA-seq pools and RNA-seq (TanyAdult_S1, TanyJuven_S2) as `trimmed_<sample>.tar` (uncompressed) and `.tar.gz`. Guppy 5.0.11 basecalls and raw Illumina reads were deleted from `~/_data/{ont,midgenomes/_sra}/` on 2026-10-23; the owner keeps copies elsewhere |
| Cluster: `/staging/l/lnell/` (was `/staging/lnell/` until 2026-10-13) | Cluster inputs and outputs. Older outputs also exist on ResearchDrive (`pelican://chtc.wisc.edu/researchdrive/arives/CHTC/`); jobs no longer write there (2026-10-14). Since 2026-10-30, `dna/trimmed/` holds only the fastp reports (`trimmed_<pool>.tar.gz`, read by `03-pool-qc`); the trimmed reads (`.tar`) were deleted for space and are on lan_images4 (`ill/dna/trimmed/`). Upload them again before rerunning `02-bwa` or `00-references/01-mitogenome` |
| Cluster: `/home/lnell/tany_time` | Clone of this repo on the submit server; jobs are submitted from it. Edit locally, push, then `git pull` there (`_bash/**/logs/` is git-ignored) |

## Conventions

- R: tidyverse, native pipe `|>`, `\(x)` lambdas, 4-space indent,
  `#'` comment blocks for prose, `# ---- ` section headers with `=*` rules.
- Derived `.rds` saved with `compress = "gz"`.
- Derived datasets are written unconditionally (since 2026-10-13), so they
  match the code after a rerun. Only slow caches are written if absent
  (`_data/bed-summs.rds` in `_R/sample-qc.R`; `allele-freqs-time-crash.rds`
  in legacy `_R/crash-test.R`); delete them to rebuild.
- Contigs are ordered as `factor(contig, levels = paste0("CONTIG", 1:45))`.
- Pool names encode site, generation, and two-digit year, e.g.
  `SN-2-07_S31` = SN, generation 2 (`b`), 2007. Years < 70 are 2000s.
- `.Rprofile` only sets the macOS quartz device size.
- `_R/pool-qc.R` needs the `rrcov` package (installed locally 2026-10-10).
  `poolfstat` is only in the container, not installed locally.
- No trailing whitespace in any file written or edited (including script
  outputs; check with `grep -n "[[:space:]]$"`).
- After every prompt, check whether `CLAUDE.md`, `PROJECT_INDEX.md`, and
  `TODO.md` need updating, unless the owner says not to.
