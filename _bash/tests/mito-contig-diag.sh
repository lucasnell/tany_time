#!/bin/bash

#' Local check (not a cluster job): what is in the Flye round-1 contig made
#' by `mito-flye-round1.sh` (with `--genome-size 16k --asm-coverage 150`)?
#'
#' Reads `_data/mito-flye-cache/with_opts_assembly.fasta` and
#' `reads1.fastq` (run `mito-flye-round1.sh` first). Reports, with
#' base-level alignments (minimap2 -c): where the two T. gracilentus COI
#' barcodes align; self-alignments off the diagonal (repeated segments);
#' and mean read depth in 1-kb bins after mapping the recruited reads back.
#' Versions are read from `_apptainer/mito-env.yml`. Requires mamba.
#' Takes ~5 min. Run from the repo root. Set ENV_DIR to reuse an
#' installed environment.
#'
#' Output: `_data/mito-contig-diag.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-contig-diag.txt
YML=${ROOT}/_bash/_apptainer/mito-env.yml
COI=${ROOT}/_bash/references/tgraci-coi.fasta
CACHE=${ROOT}/_data/mito-flye-cache
ASM=${CACHE}/with_opts_assembly.fasta
THREADS=$(sysctl -n hw.ncpu 2> /dev/null || nproc)
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/env}
cd ${WORK}

if [ ! -s ${ASM} ] || [ ! -s ${CACHE}/reads1.fastq ]; then
    echo "run _bash/tests/mito-flye-round1.sh first" 1>&2
    exit 1
fi

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

minimap2 -a -x map-ont --secondary=no -t ${THREADS} ${ASM} \
    ${CACHE}/reads1.fastq 2> /dev/null | samtools sort -o r1.bam - 2> /dev/null
samtools index r1.bam

{
    echo "minimap2 $(minimap2 --version), samtools $(pin samtools);" \
        "checked $(date +%Y-%m-%d)"
    echo "contigs (name, length):"
    awk '/^>/ {if (n) print "  " n, l; n = substr($1, 2); l = 0; next}
         {l += length($0)} END {print "  " n, l}' ${ASM}
    echo
    echo "COI barcode hits (barcode, barcode bp, contig, start, end, strand," \
        "aligned bp, identity):"
    minimap2 -c -x asm20 --secondary=yes -N 5 ${ASM} ${COI} 2> /dev/null \
        | awk '{printf "  %s %d %s %d %d %s %d %.3f\n", $1, $2, $6, $8, $9,
                $5, $4 - $3, $10 / $11}'
    echo
    echo "self-alignments >= 500 bp off the diagonal (start-end ~" \
        "matching start-end, strand, aligned bp, identity):"
    minimap2 -c -x asm20 -DP --dual=no ${ASM} ${ASM} 2> /dev/null \
        | awk '$1 == $6 && $11 >= 500 {
                printf "  %s %d-%d ~ %d-%d %s %d %.3f\n", $1, $3, $4, $8, $9,
                $5, $11, $10 / $11}'
    echo
    echo "recruited reads mapped: $(samtools view -c -F 0x904 r1.bam) of" \
        "$(awk 'END {print NR / 4}' ${CACHE}/reads1.fastq)"
    echo "mean read depth in 1-kb bins (contig, bin start kb, depth):"
    samtools depth -a r1.bam \
        | awk '{b = int(($2 - 1) / 1000); k = $1 SUBSEP b; s[k] += $3; n[k]++;
                if (!(k in seen)) {seen[k] = 1; o[++m] = k}}
               END {for (i = 1; i <= m; i++) {split(o[i], x, SUBSEP);
                    printf "  %s %d %d\n", x[1], x[2], s[o[i]] / n[o[i]]}}'
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
