#!/bin/bash

#' Local check (not a cluster job): is the whole 29.5-kb round-1 contig
#' mitochondrial? Illumina reads were not recruited with the COI bait, so
#' DNA-seq depth along the contig is unbiased by recruitment: mitochondrial
#' sequence should have depth far above nuclear, and nuclear sequence
#' (e.g., NUMT flanks) depth near the nuclear expectation.
#'
#' Maps the first N_PAIRS read pairs of two Illumina libraries to the
#' cached contig: MyKS-19-B_S18 (DNA, 2019 Myvatn pool) and TanyAdult_S1
#' (RNA-seq of adult T. gracilentus from Myvatn; depth reflects
#' transcription, so it says nothing about untranscribed sequence such as
#' the control region, but its identity to transcribed genes and COI is
#' informative). Reads are the fastp-trimmed FASTQs, streamed from the
#' uncompressed tars on the external drive (`/Volumes/lan_images4/`; the
#' 2026-10-14 output used the raw reads, since deleted locally). Contig file:
#' `_data/mito-flye-cache/with_opts_assembly.fasta` (from
#' `mito-flye-round1.sh`) with minimap2 -x sr, and reports mean depth
#' (all mapping qualities) in 1-kb bins, plus, for DNA, the nuclear
#' expectation (library bases / T. gracilentus assembly length); identity
#' (1 - NM / aligned bp) of MAPQ >= 20 reads within the gene region; and
#' identity of reads aligned to each T. gracilentus COI barcode
#' (`tgraci-coi.fasta`; the final check in `01-mitogenome` requires >= 98%
#' for both). Contig regions
#' (`mito-contig-vs-ref.sh`, `mito-contig-diag.sh`): ~3.6-19.5 kb gene
#' region; 0-3.6 and 19.5-29.5 kb control-region side, with a ~2.8-kb
#' tandem repeat in ~21.7-29.5 kb.
#' Versions are read from `_apptainer/mito-env.yml`. Requires mamba. Run
#' from the repo root. Set ENV_DIR to reuse an environment.
#'
#' Output: `_data/mito-illumina-depth.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-illumina-depth.txt
YML=${ROOT}/_bash/_apptainer/mito-env.yml
ASM=${ROOT}/_data/mito-flye-cache/with_opts_assembly.fasta
COI=${ROOT}/_bash/references/tgraci-coi.fasta
DRIVE=/Volumes/lan_images4
ASM_TG=${HOME}/_data/midgenomes/assemblies/Tgraci_assembly.fasta.gz
LIBS=(MyKS-19-B_S18:DNA:${DRIVE}/ill/dna/trimmed
      TanyAdult_S1:RNA:${DRIVE}/rna/trimmed)
N_PAIRS=4000000
THREADS=$(sysctl -n hw.ncpu 2> /dev/null || nproc)
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/env}
cd ${WORK}

[ -s ${ASM} ] || { echo "run _bash/tests/mito-flye-round1.sh first" 1>&2; exit 1; }
for x in ${LIBS[@]}; do
    [ -s ${x##*:}/trimmed_${x%%:*}.tar ] \
        || { echo "missing ${x##*:}/trimmed_${x%%:*}.tar (drive mounted?)" 1>&2; exit 1; }
done

pin () { grep -E "^\s+- $1=" ${YML} | sed 's/.*=//'; }

if [ ! -d ${ENV_DIR} ]; then
    mamba create -y -q -p ${ENV_DIR} -c conda-forge -c bioconda \
        minimap2=$(pin minimap2) samtools=$(pin samtools) > /dev/null 2>&1 \
        || CONDA_SUBDIR=osx-64 mamba create -y -q -p ${ENV_DIR} \
            -c conda-forge -c bioconda minimap2=$(pin minimap2) \
            samtools=$(pin samtools) > /dev/null 2>&1 \
        || { echo "could not install minimap2 and samtools" 1>&2; exit 1; }
fi
export PATH=${ENV_DIR}/bin:${PATH}

GL=$(gunzip -c ${ASM_TG} | awk '!/^>/ {n += length($0)} END {print n}')

{
    echo "minimap2 $(minimap2 --version), samtools $(pin samtools);" \
        "checked $(date +%Y-%m-%d)"
    echo "first ${N_PAIRS} read pairs per library; T. gracilentus assembly" \
        "${GL} bp"
    for x in ${LIBS[@]}; do
        lib=${x%%:*}
        type=$(echo ${x} | cut -d : -f 2)
        TAR=${x##*:}/trimmed_${lib}.tar
        for k in 1 2; do
            tar -xOf ${TAR} trimmed_${lib}_L002_R${k}_001.fastq.gz 2> /dev/null \
                | gunzip -c 2> /dev/null | head -n $(( N_PAIRS * 4 )) > r${k}.fq
        done
        B=$(awk 'NR % 4 == 2 {b += length($0)} END {printf "%.0f", b}' r1.fq r2.fq)
        minimap2 -a -x sr -t ${THREADS} ${ASM} r1.fq r2.fq 2> /dev/null \
            | samtools sort -o ${lib}.bam - 2> /dev/null
        samtools index ${lib}.bam
        echo
        if [ "${type}" == "DNA" ]; then
            EXP="nuclear expectation $(awk -v b=${B} -v g=${GL} \
                'BEGIN {printf "%.1f", b / g}')x"
        else
            EXP="RNA-seq, depth = transcription"
        fi
        echo "== ${lib} (${type}): ${B} bp read; ${EXP};" \
            "reads mapped $(samtools view -c -F 0x904 ${lib}.bam)," \
            "with MAPQ >= 20 $(samtools view -c -F 0x904 -q 20 ${lib}.bam)"
        echo "mean depth in 1-kb bins (bin start kb: depth):"
        samtools depth -a -Q 0 ${lib}.bam \
            | awk '{b = int(($2 - 1) / 1000); s[b] += $3; n[b]++}
                   END {for (i = 0; i in s; i++)
                            printf "%s%2d: %6d", (i % 5 ? "  " : (i ? "\n  " : "  ")),
                                i, s[i] / n[i]; printf "\n"}'
        minimap2 -c -x sr -t ${THREADS} ${ASM} r1.fq r2.fq 2> /dev/null \
            | awk '$12 >= 20 && $8 >= 4000 && $9 <= 19000 && $11 >= 100 {
                   for (i = 13; i <= NF; i++) if ($i ~ /^NM:i:/) nm = substr($i, 6)
                   s += nm; b += $11; n++}
                   END {printf "gene region (4-19 kb), MAPQ >= 20: %d reads," \
                        " identity %.4f\n", n, 1 - s / b}'
        echo "COI barcodes (barcode: reads, identity):"
        minimap2 -c -x sr -t ${THREADS} ${COI} r1.fq r2.fq 2> /dev/null \
            | awk '$11 >= 100 {for (i = 13; i <= NF; i++)
                                   if ($i ~ /^NM:i:/) nm = substr($i, 6)
                               s[$6] += nm; b[$6] += $11; n[$6]++}
                   END {for (t in n) printf "  %s: %d, %.4f\n", t, n[t],
                        1 - s[t] / b[t]}' | sort
        rm r1.fq r2.fq ${lib}.bam*
    done
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
