#!/bin/bash

#' Local check (not a cluster job): are nanopore reads longer than a
#' typical (~16-kb) mitogenome mitochondrial, and what do they contain?
#'
#' Uses the cached round-1 reads and contig from `mito-flye-round1.sh`
#' (`_data/mito-flye-cache/`). Positions refer to that 29.5-kb contig,
#' where 3600-19500 is the gene-containing segment (rRNA side to ND2 side;
#' see `mito-contig-vs-ref.sh`) and 9674-15798 is the part matched at
#' 99.9% by the high-coverage round-2 contig (`mito-round1-options.sh`).
#' Reports:
#'   1. read-length tails of the whole library (pass reads >= 1 kb, from
#'      the sequencing summary) and of the recruited reads
#'   2. nuclear coverage (pass bases / T. gracilentus assembly length)
#'   3. identity to 9674-15798 (1 - NM / block length; alignments >= 2 kb)
#'      by read length, and for reads > 17 kb with basecall q-scores
#'   4. for reads > 17 kb, alignments (>= 1 kb) to 3600-19500, classified:
#'      dimer (two same-strand copies), fold-back (two opposite-strand
#'      copies), internal (one copy with >= 2 kb unaligned on both sides),
#'      end (one copy with < 2 kb unaligned on one side)
#'   5. reads > 10 kb containing both the rRNA-side (3700-6500) and
#'      ND2-side (17400-19500) segments, and the read bp between them
#' Versions are read from `_apptainer/mito-env.yml`. Requires mamba.
#' Run from the repo root. Set ENV_DIR to reuse an environment.
#'
#' Output: `_data/mito-long-reads.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-long-reads.txt
YML=${ROOT}/_bash/_apptainer/mito-env.yml
CACHE=${ROOT}/_data/mito-flye-cache
SUMM=${HOME}/_data/ont/basecalls_guppy-5.0.11--sequencing_summary.txt.gz
ASM_TG=${HOME}/_data/midgenomes/assemblies/Tgraci_assembly.fasta.gz
THREADS=$(sysctl -n hw.ncpu 2> /dev/null || nproc)
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/env}
cd ${WORK}

[ -s ${CACHE}/reads1.fastq ] && [ -s ${CACHE}/with_opts_assembly.fasta ] \
    || { echo "run _bash/tests/mito-flye-round1.sh first" 1>&2; exit 1; }

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

R=${CACHE}/reads1.fastq
A=${CACHE}/with_opts_assembly.fasta
C=$(head -1 ${A} | sed 's/^>//; s/ .*//')
samtools faidx ${A} ${C}:9674-15798 > core.fa
samtools faidx ${A} ${C}:3600-19500 > genes.fa
samtools faidx ${A} ${C}:3700-6500 | sed '1s/.*/>rRNA_side/' > anchors.fa
samtools faidx ${A} ${C}:17400-19500 | sed '1s/.*/>ND2_side/' >> anchors.fa
gunzip -c ${SUMM} | awk -F'\t' 'NR > 1 {print $2, $10, $14, $15}' > summ.txt
awk 'NR % 4 == 1 {h = substr($1, 2)} NR % 4 == 2 {s = $0}
     NR % 4 == 0 && length(s) > 17000 {print "@" h; print s; print "+"; print}' \
    ${R} > long17.fq
awk 'NR % 4 == 1 {h = substr($1, 2)} NR % 4 == 2 {s = $0}
     NR % 4 == 0 && length(s) > 10000 {print "@" h; print s; print "+"; print}' \
    ${R} > long10.fq

tails () {  # stdin: lengths; $1: label
    awk -v lab="$1" '{l = $1; n++; b += l
        for (i = 1; i <= 5; i++) if (l > t[i]) {c[i]++; cb[i] += l}}
        BEGIN {split("10000 15000 17000 20000 25000", t, " ")}
        END {printf "  %-24s n %8d |", lab, n
             for (i = 1; i <= 5; i++) printf " %5.2f/%5.2f", 100 * c[i] / n,
                 100 * cb[i] / b; printf "\n"}'
}
nm_id () {  # PAF stdin -> read, read length, identity
    awk '{for (i = 13; i <= NF; i++) if ($i ~ /^NM:i:/) nm = substr($i, 6)
          print $1, $2, 1 - nm / $11}'
}

{
    echo "minimap2 $(minimap2 --version), samtools $(pin samtools);" \
        "checked $(date +%Y-%m-%d)"
    echo
    echo "1. % of reads / % of bases in reads longer than" \
        "10 | 15 | 17 | 20 | 25 kb"
    awk '$2 == "TRUE" && $3 >= 1000 {print $3}' summ.txt \
        | tails "library (pass, >= 1 kb)"
    awk 'NR % 4 == 2 {print length($0)}' ${R} | tails "recruited (COI bait)"
    echo
    PB=$(awk '$2 == "TRUE" {b += $3} END {printf "%.0f", b}' summ.txt)
    GL=$(gunzip -c ${ASM_TG} | awk '!/^>/ {n += length($0)} END {print n}')
    echo "2. pass bases ${PB}; assembly ${GL} bp; nuclear coverage" \
        "$(awk -v a=${PB} -v b=${GL} 'BEGIN {printf "%.0f", a / b}')x"
    echo
    echo "3. identity to ${C}:9674-15798 by read length (n, mean, median," \
        "10th and 90th percentiles):"
    minimap2 -c -x map-ont --secondary=no -t ${THREADS} core.fa ${R} \
        2> /dev/null | awk '$11 >= 2000' | nm_id > id.txt
    awk '{g = ($2 > 17000) ? "c > 17 kb" : ($2 >= 5000 && $2 <= 15000) ? \
          "b 5-15 kb" : ($2 < 5000) ? "a < 5 kb" : "x 15-17 kb"; print g "\t" $3}' \
        id.txt | sort -t$'\t' -k1,1 -k2,2n \
        | awk -F'\t' '{v[$1] = v[$1] " " $2; n[$1]++; s[$1] += $2}
            END {for (g in n) {split(substr(v[g], 2), a, " ")
                 printf "  %-10s n %4d %.4f %.4f %.4f %.4f\n", substr(g, 3),
                 n[g], s[g] / n[g], a[int(n[g] / 2) + 1], a[int(n[g] * 0.1) + 1],
                 a[int(n[g] * 0.9) + 1]}}' | sort
    echo "  mean q-score: 5-15 kb $(awk '$2 >= 5000 && $2 <= 15000 {print $1}' \
        id.txt | sort > m.n; sort -k1,1 summ.txt | join m.n - \
        | awk '{s += $4; n++} END {printf "%.2f", s / n}')," \
        "> 17 kb $(awk '$2 > 17000 {print $1}' id.txt | sort > l.n; \
        sort -k1,1 summ.txt | join l.n - | awk '{s += $4; n++} END {printf "%.2f", s / n}')"
    echo "  reads > 17 kb, identity(q-score), sorted:"
    awk '$2 > 17000 {print $1, $3}' id.txt | sort -k1,1 \
        | join - <(sort -k1,1 summ.txt) | sort -k2,2n \
        | awk '{printf " %.3f(%.1f)", $2, $5} END {print ""}' | fold -s -w 76 \
        | sed 's/^/  /; s/[[:space:]]*$//'
    echo
    echo "4. reads > 17 kb: alignments >= 1 kb to ${C}:3600-19500"
    minimap2 -x map-ont -t ${THREADS} genes.fa long17.fq 2> /dev/null \
        | awk '$11 >= 1000' | sort -k1,1 -k3,3n > g.paf
    awk 'NR % 4 == 2 {n++} END {print "  reads > 17 kb:", n}' long17.fq
    awk '{n[$1]++; L[$1] = $2; s[$1] = s[$1] $5; if (!($1 in lo)) lo[$1] = $3
          hi[$1] = $4}
         END {for (r in n) {
                if (n[r] >= 2) k = (s[r] ~ /^(\++|-+)$/) ? "dimer" : "fold-back"
                else k = (lo[r] >= 2000 && L[r] - hi[r] >= 2000) ? "internal" : "end"
                c[k]++}
              for (k in c) printf "  %s: %d\n", k, c[k]}' g.paf | sort
    echo "  with no alignment: $(( $(awk 'NR % 4 == 2' long17.fq | wc -l) - \
        $(cut -f 1 g.paf | sort -u | wc -l) ))"
    echo
    echo "5. reads > 10 kb with both rRNA-side and ND2-side segments" \
        "(read, length, bp between, order on read):"
    minimap2 -x map-ont --secondary=no -t ${THREADS} anchors.fa long10.fq \
        2> /dev/null | awk '$11 >= 1000 {print $1, $2, $6, $3, $4}' > an.txt
    awk '{L[$1] = $2; s[$1, $3] = $4; e[$1, $3] = $5; h[$1, $3]++}
         END {for (r in L) if (h[r, "rRNA_side"] && h[r, "ND2_side"]) {
                a1 = s[r, "rRNA_side"]; a2 = e[r, "rRNA_side"]
                b1 = s[r, "ND2_side"]; b2 = e[r, "ND2_side"]
                if (a2 < b1) printf "  %s %d %d rRNA_side..ND2_side\n", r, L[r], b1 - a2
                else printf "  %s %d %d ND2_side..rRNA_side\n", r, L[r], a1 - b2}}' \
        an.txt | sort -k3,3n
    echo "  (rRNA_side..ND2_side on the + strand spans the gene region," \
        "~10.9 kb on the contig)"
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
