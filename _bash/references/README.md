# `references/`: small reference sequences for cluster jobs

Written by `_R/get-mito-refs.R` (run from the repo root; queries NCBI).
Choices are logged in `_data/coi-barcode-selection.csv`.

| File | Contents | Used by |
|---|---|---|
| `tgraci-mito-bait.fasta` | All *T. gracilentus* mitochondrial records in NCBI (143; mostly COI barcodes) | `01-variant-calling/00a-mitogenome` (recruiting nanopore reads) |
| `tgraci-coi.fasta` | *T. gracilentus* COI barcodes MZ659759.1 and MN679642.1 | `00a-mitogenome` (validating the assembly) |
| `coi-barcodes.fasta` | One COI barcode per abundant trap taxon, named `bc_<code>__<accession>`: *Micropsectra lindrothi*, *M. atrofasciata* (Icelandic specimens), *Cricotopus sylvestris*, *Chironomus islandicus*, *Eukiefferiella minor*, *Simulium vittatum*, and 7 Fennoscandian *Procladius* species as a genus-level proxy (`bc_Procla`; *P. islandicus* has no NCBI sequences) | `01-variant-calling/00b-hologenome` (mitochondrial contamination index) |

Taxa follow Einarsson et al. 2026 (*Aquatic Ecology* 60:58, Table 1).
Barcode rule: COI 600–700 bp; prefer Icelandic specimens; then longest;
then lowest accession.
