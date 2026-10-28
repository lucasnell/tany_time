#!/bin/bash

#' Local check (not a cluster job): do the other species' assemblies in
#' `00-references/02-hologenome` include their mitogenomes?
#'
#' For each assembly (`~/_data/midgenomes/assemblies/`, the same files as on
#' staging), lists sequences whose names say mitochondrion, and sequences
#' the T. gracilentus gene region (`genes` in `tgraci-mito-features.fasta`)
#' matches over >= 3 kb in total (blastn, dust off, e-value <= 1e-10;
#' different subfamilies, so identity is low). A mitogenome shows up as one
#' ~15-20-kb, A+T-rich sequence covering most of the gene region. Needs
#' blastn and makeblastdb on PATH. Run from the repo root.
#'
#' Output: `_data/other-mito-check.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/other-mito-check.txt
ASM_DIR=${HOME}/_data/midgenomes/assemblies
FEATS=${ROOT}/_bash/references/tgraci-mito-features.fasta
OTHERS=("Cripar:Chironomus riparius" "Ctenta2023:Chironomus tentans (2023)"
        "Saterr:Smittia aterrima" "Pakamu:Propsilocerus akamusi")
WORK=$(mktemp -d)
cd ${WORK}
awk '/^>/ {p = ($1 == ">genes")} p' ${FEATS} > genes.fa

{
    echo "$(blastn -version | head -1); checked $(date +%Y-%m-%d)"
    for x in "${OTHERS[@]}"; do
        code=${x%%:*}
        gunzip -c ${ASM_DIR}/${code}_assembly.fasta.gz \
            | awk '/^>/ {print $1; next} {print}' > a.fa
        echo
        echo "== ${x#*:} (${code}): $(grep -c '^>' a.fa) sequences"
        echo "  named mitochondrion: $(gunzip -c ${ASM_DIR}/${code}_assembly.fasta.gz \
            | grep '^>' | grep -i 'mitochondri' | sed 's/^>//' | paste -sd ';' -)"
        makeblastdb -in a.fa -dbtype nucl -out db > /dev/null
        blastn -dust no -db db -query genes.fa -evalue 1e-10 -num_threads 8 \
            -outfmt "6 sseqid qstart qend length pident" > hits.txt
        # per sequence: aligned bp, gene-region bp covered, mean identity;
        # with the sequence's length and A+T fraction
        awk '{b[$1] += $4; id[$1] += $4 * $5
              for (i = $2; i <= $3; i++) cov[$1, i] = 1}
             END {for (s in b) {c = 0; for (k in cov) {split(k, z, SUBSEP)
                      if (z[1] == s) c++}
                  if (b[s] >= 3000) printf "%s %d %d %.1f\n", s, b[s], c, id[s] / b[s]}}' \
            hits.txt | sort -k3,3nr > top.txt
        if [ -s top.txt ]; then
            echo "  sequences matching the gene region over >= 3 kb" \
                "(name, length, A+T, aligned bp, gene-region bp covered of 15,901, mean identity):"
            while read s b c id; do
                read len at <<< $(awk -v s=${s} '/^>/ {p = (substr($1, 2) == s); next}
                    p {n += length($0); t += gsub(/[ATat]/, "")} END {printf "%d %.3f", n, t / n}' a.fa)
                printf "    %-16s %10d %.3f %6d %6d %5.1f%%\n" ${s} ${len} ${at} ${b} ${c} ${id}
            done < top.txt
        else
            echo "  no sequence matches the gene region over >= 3 kb"
        fi
        rm a.fa db.* hits.txt top.txt
    done
} > ${OUT}
sed -i '' 's/[[:space:]]*$//' ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
