#!/bin/bash

#' Collapse an assembled T. gracilentus mitochondrial contig to one copy of
#' the control-region repeat unit.
#'
#' Usage: mito-collapse.sh <contig.fasta> <features.fasta> <out prefix>
#'
#' The control region is > 22 kb and its repeat-unit copy number can't be
#' resolved from the reads (`_bash/tests/mito-copy-number.sh`), so the
#' mapping reference keeps what is known: the gene region, the unique
#' control-region sequence on each side of it, and one repeat unit.
#'
#' Features (`genes`, `unit`; from `tgraci-mito-features.fasta`) are aligned
#' to the first sequence of <contig.fasta> with minimap2 -c -x asm20. The
#' contig is reverse-complemented if the gene region is on the minus
#' strand. Pieces, on the contig:
#'   genes          best hit of `genes` (must cover >= 90% of it)
#'   unique_right   after the gene region, up to the nearest unit hit
#'                  (>= 300 bp aligned) or the contig's end
#'   unit           the unit hit with the most matching bases among those
#'                  covering >= 90% of `unit`
#'   unique_left    before the gene region, back to the nearest unit hit or
#'                  the contig's start
#' Output, in molecule order (gene region, then around the control region
#' back to the gene region's start):
#'   <prefix>.fasta       genes + unique_right + unit + unique_left, one
#'                        sequence named by the prefix's file name
#'   <prefix>_layout.tsv  piece, start and end on the contig and on the
#'                        output (1-based, inclusive), length, and whether a
#'                        side piece reached the contig's end (truncated)
#' Exits 1 without output if the gene region or a full unit is missing.
#' Needs minimap2 and samtools on PATH.
#'

CONTIG=$1
FEATS=$2
PREFIX=$3
if [ -z "${PREFIX}" ] || [ ! -s "${CONTIG}" ] || [ ! -s "${FEATS}" ]; then
    echo "Usage: $0 <contig.fasta> <features.fasta> <out prefix>" 1>&2
    exit 2
fi
NAME=$(basename ${PREFIX})
TMP=$(mktemp -d)

awk '/^>/ && n++ {exit} {print}' ${CONTIG} > ${TMP}/c0.fa
C=$(head -1 ${TMP}/c0.fa | sed 's/^>//; s/ .*//')
samtools faidx ${TMP}/c0.fa
samtools faidx ${FEATS}
samtools faidx ${FEATS} genes > ${TMP}/genes.fa
samtools faidx ${FEATS} unit > ${TMP}/unit.fa

STRAND=$(minimap2 -c -x asm20 ${TMP}/c0.fa ${TMP}/genes.fa 2> /dev/null \
    | sort -k10,10nr | awk 'NR == 1 {print $5}')
if [ "${STRAND}" == "-" ]; then
    samtools faidx -i ${TMP}/c0.fa ${C} | sed "1s/.*/>${C}/" > ${TMP}/c.fa
else
    cp ${TMP}/c0.fa ${TMP}/c.fa
fi
samtools faidx ${TMP}/c.fa
L=$(cut -f 2 ${TMP}/c.fa.fai)

# gene region: best hit covering >= 90%
read GS GE <<< $(minimap2 -c -x asm20 ${TMP}/c.fa ${TMP}/genes.fa 2> /dev/null \
    | awk '$4 - $3 >= 0.9 * $2' | sort -k10,10nr \
    | awk 'NR == 1 {print $8 + 1, $9}')
if [ -z "${GS}" ]; then
    echo "mito-collapse: no gene region covering >= 90% in ${CONTIG}" 1>&2
    rm -r ${TMP}; exit 1
fi
# unit hits: start, end, matches, covers >= 90% (1/0)
minimap2 -c -x asm20 --secondary=yes -N 50 -p 0.3 ${TMP}/c.fa ${TMP}/unit.fa \
    2> /dev/null \
    | awk '$9 - $8 >= 300 {print $8 + 1, $9, $10, ($4 - $3 >= 0.9 * $2)}' \
    | sort -k1,1n > ${TMP}/units.txt
read US UE <<< $(awk '$4 == 1' ${TMP}/units.txt | sort -k3,3nr \
    | awk 'NR == 1 {print $1, $2}')
if [ -z "${US}" ]; then
    echo "mito-collapse: no unit copy covering >= 90% in ${CONTIG}" 1>&2
    rm -r ${TMP}; exit 1
fi
# unique sides: from the gene region to the nearest unit hit
RE=$(awk -v ge=${GE} '$1 > ge {print $1 - 1; exit}' ${TMP}/units.txt)
RT=no; [ -z "${RE}" ] && { RE=${L}; RT=yes; }
LS=$(awk -v gs=${GS} '$2 < gs {s = $2 + 1} END {if (s) print s}' ${TMP}/units.txt)
LT=no; [ -z "${LS}" ] && { LS=1; LT=yes; }

{
    echo -e "piece\tcontig_start\tcontig_end\tout_start\tout_end\tlength\ttruncated"
    echo "genes ${GS} ${GE} no"
    [ ${RE} -gt ${GE} ] && echo "unique_right $((GE + 1)) ${RE} ${RT}"
    echo "unit ${US} ${UE} no"
    [ ${LS} -lt ${GS} ] && echo "unique_left ${LS} $((GS - 1)) ${LT}"
} | awk 'BEGIN {OFS = "\t"} NR == 1 {print; next}
         {n = $3 - $2 + 1; print $1, $2, $3, o + 1, o + n, n, $4; o += n}' \
    > ${PREFIX}_layout.tsv

echo ">${NAME}" > ${PREFIX}.fasta
awk -v c=${C} 'NR > 1 {print c ":" $2 "-" $3}' ${PREFIX}_layout.tsv \
    | xargs samtools faidx ${TMP}/c.fa | grep -v "^>" | tr -d '\n' \
    | fold -w 80 >> ${PREFIX}.fasta
echo >> ${PREFIX}.fasta
sed -i.bak '/^$/d' ${PREFIX}.fasta && rm ${PREFIX}.fasta.bak

rm -r ${TMP}
exit 0
