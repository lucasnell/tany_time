# `02-allele-freqs/`

mpileup files → SNAPE-pooled allele frequencies, masked allele counts, and
diversity statistics. Follows the DEST pipeline (Kapun et al. 2021,
doi:10.1093/molbev/msab259) using its helper scripts (in the container's
`/app/bin`).

| Step | Does | Jobs | Output (`/staging/lnell/dna/`) |
|---|---|---|---|
| `01-snape` | Uppercase reference and pileup; counts, coverage, and indel positions (`Mpileup2Sync.py`, base quality ≥ 25); SNAPE-pooled per contig; convert to gSYNC (`SNAPE2SYNC.py`); mask (`MaskSYNC_snape.py`) twice (below) | one per pool (`../samples/all-pools-n-adults.txt`) | `snape/<pool>_snape.tar` |
| `02-npstat` | Convert the partially masked gSYNC to pileup; npstat per contig in 10-kb windows (π, Watterson's θ, Tajima's D) at sites in the SNP file | one per pool | `npstat/<pool>_npstat_10kb.stat.gz` |
| `03-combine-syncs` | Combine fully masked gSYNC files across QC-passing pools for time, space, and all pools | one | `snape/{snape,count}_masked_{time,space,all}.tar` |

## SNAPE settings (`01-snape`)

- `-nchr` = 2 × `n_adults`; θ = 0.0088 (from `01-variant-calling/04-ont-pi`);
  D = 0.01; informative prior; unfolded.
- Masking (`MaskSYNC_snape.py`) removes sites with coverage < 10 or above
  the 95th percentile of the contig, within repeats
  (`Tgraci_repeats_locs.gff3.gz`), or within 5 bp of an indel
  (as in Kapun et al. 2021).
- **Fully masked** (`<pool>_snape_masked.sync`): also removes sites with
  P(polymorphic) < 0.9. This differs from DESTv3, which removes only
  0.1 < P < 0.9 (i.e., it keeps confidently monomorphic sites); the local
  version follows the GonzalezLab SNP-caller benchmark (see comments in
  `MaskSYNC_snape.py`). Positions kept here form the npstat SNP file.
- **Partially masked** (`<pool>_snape_part_masked.sync`): no P(polymorphic)
  filter; input to npstat so that monomorphic sites count toward diversity.

## Combined outputs (`03-combine-syncs`)

Each tar has `<name>.sync.gz`, `<name>.names.gz` (pool order), and
`<name>_noblanks.sync.gz` (sites with data in every pool).

- `snape_masked_*`: SNAPE info per pool (`ref count:alt count:ref quality:alt
  quality:P(polymorphic):P(fixed):mean minor-allele frequency`). Used by
  `_R/tidy-freq-data.R` and `_R/make-time-series-data.R`.
- `count_masked_*`: allele counts per pool (`A:T:C:G:N:del`). The `time`
  `_noblanks` file is the BayPass input (`03-select-time/00-baypass-inputs`).

Because the inputs are the fully masked files, `_noblanks` keeps only sites
that SNAPE calls polymorphic (P ≥ 0.9) in **every** pool. A site that drifted
to fixation or near-loss in any one pool is excluded (see `TODO.md`).
