#!/bin/bash

#' Local check (not a cluster job): how do the earlier whole-genome
#' assemblies of the same nanopore reads (Flye, NECAT, NextDenovo,
#' SMARTdenovo; `/Volumes/lan_images4/assembly/`, copied from
#' `/Volumes/TanyGen5TB/assembly/` 2026-10-23) lay out the mitogenome?
#'
#' Finds contigs carrying the T. gracilentus COI barcode MZ659759.1, then
#' maps onto each: the gene region of the cached round-1 contig
#' (`_data/mito-flye-cache/with_opts_assembly.fasta`; rRNA side to ND2
#' side) and its two control-region-side pieces ("CR_left", "CR_right"),
#' located by alignment (`mito-features.sh`; on the 29.5-kb Guppy contig
#' 3600-19500, 1-3600, and 19500-end, where 1-2472 is one copy of the
#' ~2.5-kb repeat unit). Gene-region copies on one contig separated by a full
#' control region give a circumference estimate (distance between their
#' starts). Assemblies are unpolished. Versions are read from
#' `_apptainer/mito-env.yml`. Requires mamba and the external drive. Run
#' from the repo root. Set ENV_DIR to reuse an environment.
#'
#' Output: `_data/mito-old-assemblies.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-old-assemblies.txt
YML=${ROOT}/_bash/_apptainer/mito-env.yml
R1C=${ROOT}/_data/mito-flye-cache/with_opts_assembly.fasta
COI=${ROOT}/_bash/references/tgraci-coi.fasta
D=/Volumes/lan_images4/assembly
OLD=(flye/contigs_flye necat/contigs_necat nextdenovo/contigs_next
     smartdenovo/contigs_smart)
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/env}
cd ${WORK}

[ -s ${R1C} ] || { echo "run _bash/tests/mito-flye-round1.sh first" 1>&2; exit 1; }
[ -d ${D} ] || { echo "external drive not mounted" 1>&2; exit 1; }

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

. ${ROOT}/_bash/tests/mito-features.sh
locate_features ${R1C} || exit 1
samtools faidx r1c.fa $(region genes) | sed '1s/.*/>genes/' > genes.fa
samtools faidx r1c.fa $(region CR_left) | sed '1s/.*/>CR_left/' > cr.fa
samtools faidx r1c.fa $(region CR_right) | sed '1s/.*/>CR_right/' >> cr.fa
cp ${COI} coi_all.fa && samtools faidx coi_all.fa MZ659759.1 > coi.fa

{
    echo "minimap2 $(minimap2 --version); checked $(date +%Y-%m-%d)"
    echo "pieces of the round-1 contig (${L} bp): genes $(region genes)," \
        "CR_left $(region CR_left), CR_right $(region CR_right)"
    for a in ${OLD[@]}; do
        f=${D}/${a}.fasta.gz
        echo
        echo "== ${a}"
        [ -s ${f} ] || { echo "  (missing)"; continue; }
        gunzip -c ${f} > asm.fa
        CTGS=$(minimap2 -c -x asm20 --secondary=yes -N 5 asm.fa coi.fa \
            2> /dev/null | cut -f 6 | sort -u)
        [ -n "${CTGS}" ] || { echo "  no contig with MZ659759.1"; continue; }
        for ctg in ${CTGS}; do
            awk -v c=${ctg} '/^>/ {p = (substr($1, 2) == c)} p' asm.fa > ctg.fa
            echo "contig ${ctg}: $(awk '!/^>/ {n += length($0)} END {print n}' \
                ctg.fa) bp, A+T $(awk '!/^>/ {s = s $0} END {n = length(s);
                a = gsub(/[AaTt]/, "", s); printf "%.3f", a / n}' ctg.fa)"
            echo "  pieces >= 500 bp (piece start-end ~ contig start-end," \
                "strand, aligned bp, identity), by contig position:"
            minimap2 -c -x asm20 --secondary=yes -N 20 -p 0.3 ctg.fa \
                genes.fa cr.fa 2> /dev/null \
                | awk '$11 >= 500 {printf "    %-8s %5d-%5d ~ %6d-%6d %s %5d %.3f\n",
                       $1, $3, $4, $8, $9, $5, $11, $10 / $11}' | sort -k3,3n
            minimap2 -c -x asm20 --secondary=yes -N 5 ctg.fa genes.fa \
                2> /dev/null | awk '$4 - $3 >= 15000 {print $8}' | sort -n \
                | awk 'NR > 1 {printf "  gene-region starts %d apart" \
                       " (circumference if the path is right)\n", $1 - p}
                       {p = $1}'
        done
    done
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
