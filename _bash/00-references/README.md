# `00-references/`: references built once

Reference files used by later steps. One job per step. Each step checks its
own output where possible and fails without writing output if the check
fails; `../check-step.sh 00-references/<step>` summarizes any step.

| Step | Does | Input (`/staging/l/lnell/`) | Output (`/staging/l/lnell/`) |
|---|---|---|---|
| `01-mitogenome` | Assemble the *T. gracilentus* mitogenome from nanopore reads of one identified male (R9.4.1, Dorado 0.9.6 sup v3.6 from `00-basecall/`): recruit reads with NCBI *T. gracilentus* mitochondrial sequences, dropping reads that align better to the nuclear assembly (its COI-carrying contig counts as mitochondrial), Flye (two rounds), medaka (`r941_min_sup_g507`), polypolish with Illumina reads from pool SN-1-08_S24. **Self-check:** fails unless both *T. gracilentus* COI barcodes (`../references/tgraci-coi.fasta`) align over ≥ 95% of their length at ≥ 98% identity | `dna/ont/`, `dna/trimmed/` (from `01-variant-calling/01-fastp`), `Tgraci_assembly.fasta.gz`, `../references/` | `Tgraci_mitogenome.tar` |
| `02-hologenome` | Combine `Tgraci_assembly.fasta.gz` with other assemblies (`nt_<code>__`: *Chironomus riparius*, *C. tentans* (2023), *Smittia aterrima*, *Propsilocerus akamusi*, the closest available genomes to abundant Mývatn Chironomini and Orthocladiinae), the *T. gracilentus* mitogenome (`mt_Tgraci__`), and COI barcodes of trap taxa (`bc_<code>__`); bwa index, `.fai`, `.dict`, and a BED of nuclear target contigs. **Self-check:** fails unless every group is present | `Tgraci_assembly.fasta.gz`, `other-genomes/`, `Tgraci_mitogenome.tar`, `../references/coi-barcodes.fasta` | `hologenome.tar` |
| `03-silva-db` | Build a Kraken2 database from SILVA 138.2 SSU rRNA (SSURef NR99; all domains, including eukaryotes) | `silva/` (three SILVA files; URLs in the script) | `kraken2-db/silva-138.2.tar.gz` |
| `04-ont-pi` | Map nanopore reads from one individual (minimap2), call variants (bcftools), mean π in 1-kb windows (vcftools). The 2021 result (0.0088; Guppy reads) is the SNAPE prior θ in `02-allele-freqs/01-snape.sh`; `../check-step.sh` compares a new estimate with it. Optional | `dna/ont/basecalls_dorado-0.9.6_sup.fastq.gz` (was Guppy 5.0.11 hac until 2026-10-23), `Tgraci_assembly.fasta.gz` | `dna/ont/ont_variants_<basecalls>.vcf.gz`, `dna/ont/ont_variants_<basecalls>.windowed.pi.gz` (the Guppy run's are untagged) |

Order: `01-mitogenome` needs trimmed reads of its polishing pool from
`01-variant-calling/01-fastp`; `02-hologenome` needs `01-mitogenome`.
`03-silva-db` and `04-ont-pi` are independent.

Before running `02-hologenome`, copy the other assemblies
(`~/_data/midgenomes/assemblies/{Cripar,Ctenta2023,Saterr,Pakamu}_assembly.fasta.gz`)
to `/staging/l/lnell/other-genomes/`. Add species to `OTHERS` in the script.
