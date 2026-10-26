#!/bin/bash

#' Mitogenome features for the local `mito-*.sh` checks, so they don't
#' depend on one contig's coordinates.
#'
#' `_bash/references/tgraci-mito-features.fasta` holds feature sequences
#' cut from the Guppy 5.0.11 round-1 contig (`_data/mito-flye-cache-guppy/with_opts_assembly.fasta`,
#' 29,542 bp), whose coordinates the checks used before 2026-10-23:
#'   genes      3600-19500   gene region, rRNA side to ND2 side
#'                           (`mito-contig-vs-ref.sh`)
#'   unit       1-2472       one copy of the ~2.5-kb control-region repeat
#'   core       9674-15798   part matched at 99.9% by the high-coverage
#'                           round-2 contig (`mito-round1-options.sh`)
#'   rRNA_side  3700-6500    anchor at the gene region's start
#'   ND2_side   17400-19500  anchor at the gene region's end
#' Rebuild it with `bash _bash/tests/mito-features.sh make` (repo root).
#'
#' Scripts source this file and call `locate_features <contig fasta>`
#' (needs minimap2 and samtools on PATH, and ROOT set to the repo root).
#' It aligns the features to the first sequence in the file (minimap2 -c
#' -x asm20) and writes, in the working directory:
#'   r1c.fa        that sequence, reverse-complemented if the gene region
#'                 aligns to the minus strand, so it always runs rRNA side
#'                 to ND2 side
#'   features.txt  name, start, end (1-based, inclusive) on r1c.fa: each
#'                 feature's best hit (most matching bases), except unit,
#'                 the leftmost hit covering >= 90% of it; plus CR_left
#'                 (1 to the gene region's start) and CR_right (its end
#'                 to the contig's end)
#' and sets C (contig name), L (length), and GS, GE (gene region start,
#' end). Fails if the gene region has no hit covering >= 90% of it; warns
#' for the other features.
#'

if [ "$1" == "make" ]; then
    G=_data/mito-flye-cache-guppy/with_opts_assembly.fasta
    [ -s ${G} ] || { echo "missing ${G}; run from the repo root" 1>&2; exit 1; }
    python3 -I - ${G} > _bash/references/tgraci-mito-features.fasta << 'PYEOF'
import sys
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
assert len(seq) == 29542, len(seq)
feats = [("genes", 3600, 19500), ("unit", 1, 2472), ("core", 9674, 15798),
         ("rRNA_side", 3700, 6500), ("ND2_side", 17400, 19500)]
for n, s, e in feats:
    print(">%s guppy_round1_contig:%d-%d" % (n, s, e))
    x = seq[s - 1:e]
    for i in range(0, len(x), 80):
        print(x[i:i + 80])
PYEOF
    exit $?
fi

locate_features () {
    local F=${ROOT}/_bash/references/tgraci-mito-features.fasta
    awk '/^>/ && n++ {exit} {print}' $1 > ctg0.fa
    C=$(head -1 ctg0.fa | sed 's/^>//; s/ .*//')
    samtools faidx ctg0.fa
    local strand=$(samtools faidx ${F} genes \
        | minimap2 -c -x asm20 ctg0.fa - 2> /dev/null \
        | sort -k10,10nr | awk 'NR == 1 {print $5}')
    if [ "${strand}" == "-" ]; then
        samtools faidx -i ctg0.fa ${C} | sed "1s/.*/>${C}/" > r1c.fa
    else
        cp ctg0.fa r1c.fa
    fi
    rm -f ctg0.fa ctg0.fa.fai
    L=$(awk '!/^>/ {n += length($0)} END {print n}' r1c.fa)
    minimap2 -c -x asm20 --secondary=yes -N 50 -p 0.3 r1c.fa ${F} 2> /dev/null \
        | awk '$4 - $3 >= 0.9 * $2 {print $1, $8 + 1, $9, $10, $5}' \
        | sort -k1,1 -k2,2n > f.hits
    awk '$1 == "unit" && !u {print $1, $2, $3; u = 1}
         $1 != "unit" && $4 > best[$1] {best[$1] = $4; s[$1] = $2 " " $3}
         END {for (k in s) print k, s[k]}' f.hits | sort > features.txt
    grep -q "^genes " features.txt \
        || { echo "locate_features: no gene region in $1" 1>&2; return 1; }
    local n
    for n in unit core rRNA_side ND2_side; do
        grep -q "^${n} " features.txt \
            || echo "locate_features: no hit for ${n} in $1" 1>&2
    done
    GS=$(awk '$1 == "genes" {print $2}' features.txt)
    GE=$(awk '$1 == "genes" {print $3}' features.txt)
    echo "CR_left 1 ${GS}" >> features.txt
    echo "CR_right ${GE} ${L}" >> features.txt
    rm -f f.hits
}
# region <name>: "contig:start-end" of a feature on r1c.fa
region () { awk -v n=$1 -v c=${C} '$1 == n {print c ":" $2 "-" $3}' features.txt; }
