#!/bin/bash

#' Local check (not a cluster job): the CONTIG37 coverage spike in pools
#' with many mitochondrial reads.
#'
#' 1. blastn of CONTIG37:2,739,300-2,739,650 against the collapsed
#'    mitogenome (`Tgraci_mitogenome.fasta`) and the uncollapsed contig
#'    (`Tgraci_mito_contig.fasta`), both from `_data/Tgraci_mitogenome.tar`.
#' 2. Share of each pool's nuclear read bases (coverage x length, from
#'    `04-bam-cov`) in CONTIG37:2,739,300-2,739,700.
#'
#' Needs blastn and the coverage files on `/Volumes/lan_images4`. Run from
#' the repo root.
#'
#' Output: `_data/contig37-spot.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/contig37-spot.txt
COV_DIR=/Volumes/lan_images4/tany_time/bam-cov
POOLS=(SN-1-96_S60 SN-2-96_S36 H-2-06_S44 SN-2-98_S66 SN-1-13_S56
       KS-1-15-B_S12 SN-2-97_S37 Lei-19_S13)
WORK=$(mktemp -d)
cd ${WORK}

tar -xf ${ROOT}/_data/Tgraci_mitogenome.tar \
    Tgraci_mitogenome/Tgraci_mitogenome.fasta \
    Tgraci_mitogenome/Tgraci_mito_contig.fasta
cat Tgraci_mitogenome/Tgraci_mitogenome.fasta \
    Tgraci_mitogenome/Tgraci_mito_contig.fasta > mito.fa
{
    echo ">spot"
    gunzip -c ${HOME}/_data/midgenomes/assemblies/Tgraci_assembly.fasta.gz \
        | awk '/^>/ {p = ($1 == ">CONTIG37"); next} p' | tr -d '\n' \
        | cut -c2739300-2739650
    echo
} > spot.fa

{
    echo "$(blastn -version | head -1); checked $(date +%Y-%m-%d)"
    echo
    echo "== blastn: CONTIG37:2,739,300-2,739,650 (query 1 = 2,739,300)" \
        "vs mitogenome (e-value <= 1e-10)"
    echo "subject qstart qend sstart send length pident evalue"
    blastn -query spot.fa -subject mito.fa -dust no -word_size 11 \
        -evalue 1e-10 -outfmt "6 sseqid qstart qend sstart send length pident evalue" \
        | sort -k8,8g | tr '\t' ' '
    echo
    echo "== Share of nuclear read bases in CONTIG37:2,739,300-2,739,700," \
        "and maximum depth there"
    for p in ${POOLS[@]}; do
        gunzip -c ${COV_DIR}/${p}_bwa_cov.bed.gz \
            | awk -v p=${p} '{b = ($3 - $2) * $4; t += b}
                $1 == "CONTIG37" && $3 > 2739300 && $2 < 2739700 {
                    s += b; if ($4 > m) m = $4}
                END {printf "%-14s %.4f %d\n", p, s / t, m}'
    done
} > ${OUT}
sed -i '' 's/[[:space:]]*$//' ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
