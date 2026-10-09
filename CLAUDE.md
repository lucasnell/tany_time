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
  sequence; no further pools can be added.

See `PROJECT_INDEX.md` for the file map and `TODO.md` for open items.

## Layout

- `_bash/` — cluster pipeline (UW–Madison CHTC, HTCondor). Each step is a
  `.sh` + `.sub` pair, numbered in run order. Runs inside the Apptainer image
  built from `_bash/_apptainer/tany_time.def` (`/staging/lnell/tany_time.sif`).
- `_R/` — current (2025) local R analysis. Run from the repo root.
- `_data/` — small inputs and derived datasets. `*.rds`, `*.gz`, and `*.pdf`
  are git-ignored, so derived `.rds` files exist only locally.
- `_plots/` — figures (git-ignored PDFs).
- `_refs/` — reference documents (BayPass v2.31 manual; git-ignored PDF).
- `_legacy/` — older analyses (mostly 2020–2022) formerly at the repo root (previous
  MAPGD/SNAPE pipeline, spatial FST, admixture, simulations). They read from
  old paths (`~/Box Sync`, `~/Stanford_Drive`, `~/Cornell_Drive`,
  `full-DNA-info.csv`), and their relative paths (`sourceCpp("helpers.cpp")`,
  `_plots/`, `.Rprofile`) assumed the repo root, so they will not run as-is.
  Do not edit unless asked. `_R/crash-test.R` is also legacy (superseded by
  BayPass).

## Data outside the repo

Scripts read large files by absolute path; none are in git.

| Path | Contents |
|---|---|
| `~/_data/_tany-time/snape/{snape,count}_masked_{time,space,all}/` | Combined SNAPE frequency and read-count sync files |
| `~/_data/_baypass/{baypass-inputs,init,beta}/` | BayPass inputs and outputs (local copy holds only `beta0` runs 1–10) |
| `~/_data/midgenomes/` | Assembly, GO terms, OrthoFinder HOG→GO tables |
| `~/_data/Tgraci_genes.gff3` | Gene annotation |
| `~/Desktop/Tgraci_proteins.faa`, `~/Desktop/Tgraci-funct-Dmel.csv` | Protein FASTA and BLASTP output used in `baypass-analysis.R` |
| `~/_data/_tany-time/trimmed/trimmed_reports/` | fastp/FastQC per-sample reports and `multiqc_data/` (`sample-qc.R` still points to the old `~/_data/trimmed_reports`) |
| `/Volumes/Lucas512/...` | BAM coverage `.bed.gz` files (external drive) |
| Cluster: `/staging/lnell/`, `pelican://chtc.wisc.edu/researchdrive/arives/CHTC/` | Cluster inputs and outputs |

## Conventions

- R: tidyverse, native pipe `|>`, `\(x)` lambdas, 4-space indent,
  `#'` comment blocks for prose, `# ---- ` section headers with `=*` rules.
- Derived `.rds` saved with `compress = "gz"`.
- Many writes are guarded by `if (!file.exists(...))` (or an `OVERWRITE`
  flag), so a saved file can be older than the code that made it. Check
  before trusting one; `_R/baypass-check-candidates.R` does this for
  `candi_df.rds`.
- Contigs are ordered as `factor(contig, levels = paste0("CONTIG", 1:45))`.
- Pool names encode site, generation, and two-digit year, e.g.
  `SN-2-07_S31` = SN, generation 2 (`b`), 2007. Years < 70 are 2000s.
- `.Rprofile` only sets the macOS quartz device size.
- After every prompt, check whether `CLAUDE.md`, `PROJECT_INDEX.md`, and
  `TODO.md` need updating, unless the owner says not to.
