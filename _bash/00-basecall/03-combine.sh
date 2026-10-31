#!/bin/bash

#' Combine the Dorado chunk outputs (`02-dorado`) into one gzipped FASTQ of
#' reads with mean quality >= MIN_QS (10, MinKNOW's default pass threshold
#' for sup) and one summary table of all reads. Fails if any chunk listed
#' in `manifest.tsv` lacks a BAM or summary.
#'
#' Output (`/staging/l/lnell/dna/ont/`):
#' `basecalls_dorado-0.9.6_sup.fastq.gz`,
#' `basecalls_dorado-0.9.6_sup--summary.tsv.gz`.
#'


. /app/.bashrc
. bin/helpers.sh || exit 1
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export CHUNK_DIR="/staging/l/lnell/dna/ont/fast5_chunks"
export PARENT_DIR_IN="/staging/l/lnell/dna/ont/dorado"
export PARENT_DIR_OUT="/staging/l/lnell/dna/ont"

export MIN_QS=10
export OUT_FQ=basecalls_dorado-0.9.6_sup.fastq.gz
export OUT_SUMM=basecalls_dorado-0.9.6_sup--summary.tsv.gz

require_files ${CHUNK_DIR}/manifest.tsv
CHUNKS=($(cut -f 1 ${CHUNK_DIR}/manifest.tsv | sed 's/\.tar$//' | sort -u))
for c in ${CHUNKS[@]}; do
    require_files ${PARENT_DIR_IN}/${c}.bam ${PARENT_DIR_IN}/${c}_summary.tsv
done
echo "chunks: ${#CHUNKS[@]}"

: > ${OUT_FQ}
for c in ${CHUNKS[@]}; do
    samtools fastq -@ ${THREADS} -T qs ${PARENT_DIR_IN}/${c}.bam 2> /dev/null \
        | awk -v m=${MIN_QS} 'NR % 4 == 1 {keep = 0
                for (i = 2; i <= NF; i++) if ($i ~ /^qs:f:/) keep = (substr($i, 6) + 0 >= m)}
              keep' \
        | gzip >> ${OUT_FQ}
    check_exit_status "fastq ${c}" $?
done

head -n 1 ${PARENT_DIR_IN}/${CHUNKS[0]}_summary.tsv > summary.tsv
for c in ${CHUNKS[@]}; do
    tail -n +2 ${PARENT_DIR_IN}/${c}_summary.tsv >> summary.tsv
done
gzip summary.tsv && mv summary.tsv.gz ${OUT_SUMM}
check_exit_status "combine summaries" $?

# Totals: all reads from the summaries; pass reads from the FASTQ
gunzip -c ${OUT_SUMM} | awk -F'\t' -v m=${MIN_QS} '
    NR == 1 {for (i = 1; i <= NF; i++) {if ($i == "sequence_length_template") L = i
                                        if ($i == "mean_qscore_template") Q = i}; next}
    {n++; b += $L; q += $Q; if ($Q >= m) {p++; pb += $L}}
    END {printf "all reads: %d (%.0f bp, mean q %.2f); q >= %s: %d reads, %.0f bp\n",
         n, b, q / n, m, p, pb}'
echo "FASTQ reads: $(gunzip -c ${OUT_FQ} | awk 'END {print NR / 4}')"

mv ${OUT_FQ} ${OUT_SUMM} ${PARENT_DIR_OUT}/
check_exit_status "move output to staging" $?

safe_exit 0
