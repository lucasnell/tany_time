# tany_time

Pool-seq population genomics of the midge *Tanytarsus gracilentus* at Lake
Mývatn, Iceland: how population crashes affect its adaptive evolution and
neutral genetic variation.

Data are pooled whole-genome sequences of male midges from archived
emergence-trap samples (one site, 1992–2015, two generations per year) and
from 2019 samples across Mývatn and other Icelandic lakes, alongside the
long-term abundance time series from the same trap.

## Layout

| Folder | Contents |
|---|---|
| `_bash/` | Cluster pipeline (CHTC/HTCondor, Apptainer): reads → allele frequencies → BayPass. See `_bash/README.md`. |
| `_R/` | Local analyses, run from the repo root |
| `_data/` | Sample information, abundances, and small derived datasets (large files are not in git) |
| `_plots/` | Figures (not in git) |
| `_refs/` | Reference documents (not in git) |
| `_legacy/` | Earlier analyses (2020–2022), kept for reference; not runnable as-is |

## Project notes

- `CLAUDE.md`: project context, conventions, and data locations
- `PROJECT_INDEX.md`: file map and pipeline order
- `TODO.md`: open items and decisions
