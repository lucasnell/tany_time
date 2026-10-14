# `tests/`: local checks of pipeline behavior

Scripts run locally (not on the cluster) from the repo root. Each writes
its result to `_data/`.

| Script | Checks | Output |
|---|---|---|
| `npstat-maxcov.sh` | npstat output is identical when `-maxcov` is lowered to the maximum coverage present, and how much faster that is (`02-allele-freqs/02-npstat` relies on this) | `_data/npstat-maxcov-check.txt` |
| `v2-tool-check.sh` | Each pipeline tool call works with the versions pinned in `_apptainer/*.yml` (run after changing pins); Python-script, π, and BayPass outputs match or are read as expected | `_data/v2-tool-check.txt` |
| `mito-flye-round1.sh` | Flye round 1 of `00-references/01-mitogenome` as written and with `--genome-size 16k --asm-coverage 150`; caches recruited reads and assembly in `_data/mito-flye-cache/` | `_data/mito-flye-round1-check.txt` |
| `mito-contig-diag.sh` | COI barcode hits, self-repeats, and read depth along the cached round-1 contig | `_data/mito-contig-diag.txt` |
| `mito-round1-options.sh` | Round-1 alternatives (trim to high-depth core, recruit >= 1 kb, Flye `--meta`) each carried through round 2 | `_data/mito-round1-options.txt` |
| `mito-contig-vs-ref.sh` | Which parts of the T. verralli and T. formosanus mitogenomes each test contig covers | `_data/mito-contig-vs-ref.txt` |
| `mito-illumina-depth.sh` | Illumina depth along the round-1 contig (MyKS-19-B_S18 DNA pool; TanyAdult_S1 RNA-seq), and identity to the gene region and to each COI barcode | `_data/mito-illumina-depth.txt` |
| `mito-old-assemblies.sh` | Layout of the mitochondrial contigs in the 2019 whole-genome assemblies (`/Volumes/TanyGen5TB/assembly/`): gene region, control-region pieces, repeat-unit copies | `_data/mito-old-assemblies.txt` |
| `mito-long-reads.sh` | Whether recruited reads > 17 kb are mitochondrial (length tails, identity by length and q-score, read structure) | `_data/mito-long-reads.txt` |
