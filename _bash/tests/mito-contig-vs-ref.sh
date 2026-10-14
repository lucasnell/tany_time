#!/bin/bash

#' Local check (not a cluster job): which parts of the Tanytarsus
#' mitogenome do the test Flye contigs cover?
#'
#' Downloads the complete mitogenomes of Tanytarsus verralli (PV994466.1)
#' and T. formosanus (NC_069313.1) with their gene coordinates (NCBI
#' E-utilities), then aligns every contig from `mito-flye-round1.sh` and
#' `mito-round1-options.sh` to each (minimap2 -c -x asm20). Reports
#' alignment blocks (contig start-end ~ reference start-end) and the
#' contig bp left unaligned. A contig covering part of the reference twice
#' is duplicated; unaligned contig sequence is not like either reference
#' (non-mitochondrial, or a control region, which varies among species).
#' Versions are read from `_apptainer/mito-env.yml`. Requires mamba and
#' internet. Run from the repo root. Set ENV_DIR to reuse an environment.
#'
#' Output: `_data/mito-contig-vs-ref.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-contig-vs-ref.txt
YML=${ROOT}/_bash/_apptainer/mito-env.yml
CACHE=${ROOT}/_data/mito-flye-cache
REFS=(PV994466.1 NC_069313.1)
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/env}
EU="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
cd ${WORK}

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

for r in ${REFS[@]}; do
    curl -s -G ${EU} --data-urlencode db=nuccore --data-urlencode id=${r} \
        --data-urlencode rettype=fasta > ${r}.fasta
    sleep 1
    curl -s -G ${EU} --data-urlencode db=nuccore --data-urlencode id=${r} \
        --data-urlencode rettype=ft > ${r}.ft
    sleep 1
    [ -s ${r}.fasta ] && [ -s ${r}.ft ] \
        || { echo "download of ${r} failed" 1>&2; exit 1; }
    cp ${r}.fasta ${CACHE}/
done

# Contigs: label = file stem, contig name prefixed
: > contigs.fasta
for f in ${CACHE}/with_opts_assembly.fasta ${CACHE}/options/*.fasta; do
    [ -s ${f} ] || continue
    s=$(basename ${f} .fasta)
    sed "s/^>\([^ ]*\).*/>${s}:\1/" ${f} >> contigs.fasta
done
samtools faidx contigs.fasta

{
    echo "minimap2 $(minimap2 --version); checked $(date +%Y-%m-%d)"
    for r in ${REFS[@]}; do
        echo
        echo "== reference $(head -1 ${r}.fasta | sed 's/^>//')"
        echo "genes (start end name; from the NCBI feature table):"
        # gene features: "start end gene" then a "gene NAME" qualifier line
        awk -F'\t' '$3 == "gene" {s = $1; e = $2; want = 1; next}
                    want && $4 == "gene" {printf "  %s %s %s\n", s, e, $5;
                                          want = 0}' ${r}.ft \
            | tr -d '<>' | tr '\n' ';' | sed 's/;  /; /g; s/;$/\n/' \
            | fold -s -w 76 | sed 's/^/  /; s/[[:space:]]*$//'
        echo "contig blocks (contig, length: contig start-end ~ ref start-end" \
            "strand identity; unaligned bp):"
        minimap2 -c -x asm20 --secondary=yes -N 5 ${r}.fasta contigs.fasta \
            2> /dev/null | sort -k1,1 -k3,3n > hits.paf
        while read c len; do
            echo "  ${c}, ${len} bp:"
            awk -v c=${c} '$1 == c {printf "    %d-%d ~ %d-%d %s %.3f\n",
                           $3, $4, $8, $9, $5, $10 / $11}' hits.paf
            awk -v c=${c} -v L=${len} '$1 == c {for (i = $3; i < $4; i++) a[i] = 1}
                END {n = 0; for (i in a) n++; printf "    unaligned: %d bp\n", L - n}' \
                hits.paf
        done < <(cut -f 1,2 contigs.fasta.fai)
    done
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
