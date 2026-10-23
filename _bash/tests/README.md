# `tests/`: local checks of pipeline behavior

Scripts run locally (not on the cluster) from the repo root. Each writes
its result to `_data/`.

Large inputs are on the external drive `/Volumes/lan_images4/` (sup
nanopore reads; trimmed Illumina tars). The Guppy 5.0.11 files read by
`dorado-r941-test.sh`, `dorado-sup-check.sh`, and `mito-long-reads.sh`
were deleted locally on 2026-10-23, so those three need them restored
to rerun. Mitogenome outputs saved before then used the Guppy reads.

The `mito-*.sh` checks locate the gene region, repeat unit, and other
features on whatever round-1 contig is cached by aligning the sequences in
`_data/mito-features.fasta` (`mito-features.sh`; cut from the Guppy
round-1 contig, now in `_data/mito-flye-cache-guppy/`), instead of fixed
coordinates.

| Script | Checks | Output |
|---|---|---|
| `npstat-maxcov.sh` | npstat output is identical when `-maxcov` is lowered to the maximum coverage present, and how much faster that is (`02-allele-freqs/02-npstat` relies on this) | `_data/npstat-maxcov-check.txt` |
| `v2-tool-check.sh` | Each pipeline tool call works with the versions pinned in `_apptainer/*.yml` (run after changing pins); Python-script, π, and BayPass outputs match or are read as expected | `_data/v2-tool-check.txt` |
| `mito-flye-round1.sh` | Flye round 1 of `00-references/01-mitogenome` as written and with `--genome-size 16k --asm-coverage 150`; caches recruited reads and assembly in `_data/mito-flye-cache/` | `_data/mito-flye-round1-check.txt` |
| `dorado-r941-test.sh` | Dorado 0.9.6 basecalls the 2019 R9.4.1 FAST5 data with the sup model; predicted quality vs Guppy 5.0.11 hac on the same reads | `_data/dorado-r941-test.txt` |
| `dorado-sup-check.sh` | Combined sup basecalls from `00-basecall/03-combine` (local copy on `/Volumes/lan_images4/ont/`): FASTQ integrity, q >= 10 reads match the summary, yield, N50, read overlap, and predicted quality vs Guppy 5.0.11 hac | `_data/dorado-sup-check.txt` |
| `mito-reassemble.sh` | All mitochondrial reads reassembled with Flye, Raven, NextDenovo, and SMARTdenovo; control-region repeat structure tested by read bridging and Illumina copy number | `_data/mito-reassemble.txt` |
| `mito-cr-span.sh` | Control-region length per molecule: reads crossing from the gene region's end to its start (checked on synthetic reads), and how far reads run past each end without crossing | `_data/mito-cr-span.txt` |
| `mito-depth-by-length.sh` | Nanopore depth by read length in the gene region, the unique control-region segments, and the repeat copies of the reassembled SMARTdenovo contig | `_data/mito-depth-by-length.txt` |
| `mito-contig-diag.sh` | COI barcode hits, self-repeats, and read depth along the cached round-1 contig | `_data/mito-contig-diag.txt` |
| `mito-round1-options.sh` | Round-1 alternatives (trim to high-depth core, recruit >= 1 kb, Flye `--meta`) each carried through round 2 | `_data/mito-round1-options.txt` |
| `mito-contig-vs-ref.sh` | Which parts of the T. verralli and T. formosanus mitogenomes each test contig covers | `_data/mito-contig-vs-ref.txt` |
| `mito-illumina-depth.sh` | Illumina depth along the round-1 contig (MyKS-19-B_S18 DNA pool; TanyAdult_S1 RNA-seq), and identity to the gene region and to each COI barcode | `_data/mito-illumina-depth.txt` |
| `mito-old-assemblies.sh` | Layout of the mitochondrial contigs in the 2019 whole-genome assemblies (`/Volumes/TanyGen5TB/assembly/`): gene region, control-region pieces, repeat-unit copies | `_data/mito-old-assemblies.txt` |
| `mito-features-check.sh` | `locate_features` (`mito-features.sh`) recovers the fixed coordinates used on the Guppy round-1 contig; feature positions on the sup round-1 contig | `_data/mito-features-check.txt` |
| `medaka-model-test.sh` | Medaka `r941_min_hac_g507` and `r941_min_sup_g507` on the sup v3.6 reads vs no polishing: errors against Illumina (MyKS-19-B, 4M trimmed pairs cached in `_data/mito-flye-cache/ill/`) in and outside the gene region | `_data/medaka-model-test.txt` |
| `mito-long-reads.sh` | Whether recruited reads > 17 kb are mitochondrial (length tails, identity by length and q-score, read structure) | `_data/mito-long-reads.txt` |
