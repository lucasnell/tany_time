# `references/`: small reference sequences for cluster jobs

Written by `_R/get-mito-refs.R` (run from the repo root; queries NCBI).
Choices are logged in `_data/coi-barcode-selection.csv`.

| File | Contents | Used by |
|---|---|---|
| `tgraci-mito-bait.fasta` | All *T. gracilentus* mitochondrial records in NCBI (143; mostly COI barcodes) | `00-references/01-mitogenome` (recruiting nanopore reads) |
| `tgraci-coi.fasta` | *T. gracilentus* COI barcodes MZ659759.1 and MN679642.1 | `00-references/01-mitogenome` (validating the assembly) |
| `tgraci-mito-features.fasta` | Not from NCBI: gene region, one control-region repeat unit, and anchors cut from the Guppy 5.0.11 round-1 contig (`_bash/tests/mito-features.sh make`) | `00-references/01-mitogenome` (`mito-collapse.sh`), `_bash/tests/mito-*.sh` |
| `other-mitogenomes.fasta` | Mitogenomes of *C. riparius* (OU815687.1) and *S. aterrima* (NC_069312.1), whose assemblies lack one; named `mt_<code>__<accession>`; from `_R/get-other-mitogenomes.R` | `00-references/02-hologenome` |
| `coi-barcodes.fasta` | One COI barcode per abundant trap taxon, named `bc_<code>__<accession>`: *Micropsectra lindrothi*, *M. atrofasciata* (Icelandic specimens), *Cricotopus sylvestris*, *Chironomus islandicus*, *Eukiefferiella minor*, *Simulium vittatum*, and 7 Fennoscandian *Procladius* species as a genus-level proxy (`bc_Procla`; *P. islandicus* has no NCBI sequences) | `00-references/02-hologenome` (mitochondrial contamination index) |

Taxa follow Einarsson et al. 2026 (*Aquatic Ecology* 60:58, Table 1).
Barcode rule: COI 600–700 bp; prefer Icelandic specimens; then longest;
then lowest accession.
