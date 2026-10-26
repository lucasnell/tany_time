#!/bin/bash

#' Local check (not a cluster job): how long is the control region, from
#' nanopore read depth and from the absence of reads crossing it?
#'
#' 1. Depth. Builds a reference of four pieces of the round-1 contig
#'    (located by `mito-features.sh`): the gene region, one repeat unit,
#'    and the unique control-region sequence on each side of the gene
#'    region (control-region sequence not covered by hits of the unit,
#'    >= 300 bp at any identity minimap2 -x asm20 reports). Maps all
#'    filtered mitochondrial reads (`mito-reassemble.sh` cache:
#'    `_data/mito-flye-cache/reassemble/mito_reads.fastq`; minimap2 -c -x
#'    map-ont, primary and supplementary alignments, so a read spanning k
#'    unit copies aligns about k times to the one unit) and sums aligned
#'    reference bases per piece. Depth = bases / piece length. Unit copies
#'    per molecule = unit depth / gene-region depth; control-region length =
#'    (unique-left + unique-right + copies x unit length). This assumes
#'    reads sample all parts of the molecule equally; the unique pieces'
#'    depth relative to the gene region tests that (1 if equal), so the
#'    estimate is also given with unit depth relative to the unique pieces.
#'    Reported by read length (< 5, 5-11, >= 11 kb) and overall.
#' 2. Crossing. A read crosses the control region when it contains the
#'    whole control region plus >= 500 bp of the gene region on each side
#'    (as `mito-cr-span.sh` detects; none observed). If reads start
#'    uniformly on a circular molecule of length G + L (gene region G,
#'    control region L), a read of length l crosses with probability
#'    max(0, l - L - 1000) / (G + L), so the expected number of crossing
#'    reads is E(L) = sum over reads of that. With none seen, lengths with
#'    E(L) > 3 are rejected at the 5% level (Poisson P(0) = exp(-E)). This
#'    gives a lower bound only: longer control regions predict fewer
#'    crossings, so no crossings fit any length above the bound.
#' 3. Depth along a full-length reference: the same reads mapped (primary
#'    alignments, all MAPQ) to the 2019 NextDenovo contig ctg000860 from its
#'    start to its second gene-region copy (gene region, then the whole
#'    ~28.7-kb control region with ~8 unit copies; `mito-old-assemblies.sh`),
#'    mean depth in 2.5-kb bins. Tests whether depth is even along the gene
#'    region and how control-region depth compares without collapsing the
#'    repeats.
#' Versions are read from `_apptainer/mito-env.yml`. Requires mamba and
#' `/Volumes/lan_images4/assembly/` (step 3). Run from the repo root. Set
#' ENV_DIR to reuse an environment; OUT and MR override the output and
#' read files.
#'
#' Output: `_data/mito-copy-number.txt`.
#'

ROOT=$(pwd)
OUT=${OUT:-${ROOT}/_data/mito-copy-number.txt}
YML=${ROOT}/_bash/_apptainer/mito-env.yml
CACHE=${ROOT}/_data/mito-flye-cache
R1C=${CACHE}/with_opts_assembly.fasta
MR=${MR:-${CACHE}/reassemble/mito_reads.fastq}
NEXT=/Volumes/lan_images4/assembly/nextdenovo/contigs_next.fasta.gz
THREADS=$(sysctl -n hw.ncpu 2> /dev/null || nproc)
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/env}
cd ${WORK}

[ -s ${R1C} ] && [ -s ${MR} ] \
    || { echo "run _bash/tests/mito-reassemble.sh first" 1>&2; exit 1; }
[ -s ${NEXT} ] || { echo "missing ${NEXT} (drive mounted?)" 1>&2; exit 1; }

pin () { grep -E "^\s+- $1=" ${YML} | sed 's/.*=//'; }

if [ ! -x ${ENV_DIR}/bin/minimap2 ]; then
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
samtools faidx r1c.fa
samtools faidx r1c.fa $(region unit) | sed '1s/.*/>unit/' > unit.fa

# Control-region positions not covered by unit hits; longest stretch on
# each side of the gene region
samtools faidx r1c.fa ${C}:1-$((GS - 1)) | sed '1s/.*/>left/' > left.fa
samtools faidx r1c.fa ${C}:$((GE + 1))-${L} | sed '1s/.*/>right/' > right.fa
cat left.fa right.fa > cr.fa
minimap2 -c -x asm20 --secondary=yes -N 50 -p 0.3 cr.fa unit.fa 2> /dev/null \
    | awk '$9 - $8 >= 300 {print $6, $8 + 1, $9}' > unit_hits.txt
for side in left right; do
    len=$(awk '!/^>/ {n += length($0)} END {print n}' ${side}.fa)
    awk -v s=${side} -v n=${len} '
        $1 == s {for (i = $2; i <= $3; i++) m[i] = 1}
        END {for (i = 1; i <= n + 1; i++) {
                 if (i <= n && !(i in m)) {if (!st) st = i}
                 else if (st) {if (i - st > best) {best = i - st; bs = st}; st = 0}}
             if (best) print s, bs, bs + best - 1}' unit_hits.txt
done > unique.txt
{
    samtools faidx r1c.fa $(region genes) | sed '1s/.*/>genes/'
    cat unit.fa
    while read side s e; do
        samtools faidx ${side}.fa ${side}:${s}-${e} | sed "1s/.*/>unique_${side}/"
    done < unique.txt
} > ref.fa
samtools faidx ref.fa

# Aligned reference bases per piece and read-length class
minimap2 -c -x map-ont --secondary=no -t ${THREADS} ref.fa ${MR} 2> /dev/null \
    | awk '{c = ($2 < 5000) ? "a" : ($2 < 11000) ? "b" : "c"
            b[c, $6] += $9 - $8; b["all", $6] += $9 - $8}
           END {for (k in b) {split(k, x, SUBSEP); print x[1], x[2], b[k]}}' \
    > bases.txt

awk 'NR % 4 == 2 {print length($0)}' ${MR} > lens.txt

# Full-length reference: ctg000860 up to its second gene-region copy
gunzip -c ${NEXT} | awk '/^>/ {p = (substr($1, 2) == "ctg000860")} p' > next.fa
samtools faidx next.fa
samtools faidx r1c.fa $(region genes) | sed '1s/.*/>genes/' > genes.fa
minimap2 -c -x asm20 --secondary=yes -N 5 next.fa genes.fa 2> /dev/null \
    | awk '$4 - $3 >= 15000 {print $8 + 1, $9}' | sort -n > next_genes.txt
NG1=$(awk 'NR == 1 {print $1}' next_genes.txt)
NGE1=$(awk 'NR == 1 {print $2}' next_genes.txt)
NG2=$(awk 'NR == 2 {print $1}' next_genes.txt)
samtools faidx next.fa ctg000860:${NG1}-$((NG2 - 1)) | sed '1s/.*/>full/' > full.fa
minimap2 -a -x map-ont --secondary=no -t ${THREADS} full.fa ${MR} 2> /dev/null \
    | samtools sort -o full.bam - 2> /dev/null
samtools index full.bam
samtools depth -a -Q 0 full.bam \
    | awk '{b = int(($2 - 1) / 2500); s[b] += $3; n[b]++}
           END {for (i = 0; i in s; i++) printf "%d %d %d\n", i * 2500 + 1, i * 2500 + n[i], s[i] / n[i]}' \
    > full_bins.txt

{
    echo "minimap2 $(minimap2 --version), samtools $(pin samtools);" \
        "checked $(date +%Y-%m-%d)"
    echo "round-1 contig ${C} (${L} bp): gene region $(region genes)," \
        "unit $(region unit); unit hits on the control-region sides" \
        "(side start end): $(awk '{printf "%s%s %d-%d", (NR > 1 ? ", " : ""), $1, $2, $3}' unit_hits.txt)"
    echo "reference pieces (name length):" \
        "$(awk '{printf " %s %d", $1, $2}' ref.fa.fai)"
    echo "reads: $(wc -l < lens.txt | tr -d ' ') (filtered mitochondrial)"
    echo
    echo "1. depth (aligned bases / piece length) and estimates"
    echo "   class: depth genes unit unique_left unique_right | unique/gene" \
        "| copies vs gene, vs unique | control region kb vs gene, vs unique"
    for c in a b c all; do
        awk -v c=${c} 'NR == FNR {len[$1] = $2; next}
            $1 == c {d[$2] = $3 / len[$2]}
            END {lab = (c == "a") ? "< 5 kb" : ((c == "b") ? "5-11 kb" : ((c == "c") ? ">= 11 kb" : "all"))
                 g = d["genes"]; u = d["unit"]
                 ul = d["unique_left"]; ur = d["unique_right"]
                 uq = (ul * len["unique_left"] + ur * len["unique_right"]) / (len["unique_left"] + len["unique_right"])
                 kg = u / g; ku = u / uq
                 fl = len["unique_left"] + len["unique_right"]
                 printf "   %-8s %7.0f %7.0f %7.0f %7.0f | %4.2f | %5.1f %5.1f | %5.1f %5.1f\n",
                     lab, g, u, ul, ur, uq / g, kg, ku,
                     (fl + kg * len["unit"]) / 1000, (fl + ku * len["unit"]) / 1000}' \
            ref.fa.fai bases.txt
    done
    echo
    echo "2. expected reads crossing the control region, E(L), for control" \
        "region L (gene region $((GE - GS + 1)) bp); none observed"
    G=$((GE - GS + 1))
    for Lc in 15000 18000 20000 21000 22000 23000 24000 25000 27000 30000; do
        awk -v L=${Lc} -v G=${G} '{x = $1 - L - 1000; if (x > 0) e += x}
            END {printf "   L %6d: E %7.2f, P(0) %.3g\n", L, e / (G + L), exp(-e / (G + L))}' lens.txt
    done
    awk -v G=${G} '{l[NR] = $1}
        END {for (L = 10000; L <= 60000; L += 100) {e = 0
                 for (i = 1; i <= NR; i++) {x = l[i] - L - 1000; if (x > 0) e += x}
                 if (e / (G + L) <= 3) {printf "   lower bound (E <= 3, 5%%): %d bp\n", L; exit}}}' lens.txt
    echo "   longest reads (bp): $(sort -nr lens.txt | head -5 | tr '\n' ' ')"
    echo
    echo "3. mean depth (primary, all MAPQ) in 2.5-kb bins along ctg000860" \
        "${NG1}-$((NG2 - 1)) (2019 NextDenovo): gene region 1-$((NGE1 - NG1 + 1))," \
        "control region after"
    awk -v ge=$((NGE1 - NG1 + 1)) '{printf "   %6d-%6d %-15s %6d\n", $1, $2,
        ($2 <= ge ? "gene region" : ($1 > ge ? "control region" : "both")), $3}' full_bins.txt
} > ${OUT}
sed -i '' 's/[[:space:]]*$//' ${OUT} 2> /dev/null || sed -i 's/[[:space:]]*$//' ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
