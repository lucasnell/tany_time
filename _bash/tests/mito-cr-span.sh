#!/bin/bash

#' Local check (not a cluster job): how long is the control region in the
#' molecules the nanopore reads came from?
#'
#' The gene region (feature `genes` of `mito-features.sh`, located on the
#' round-1 contig: rRNA side to ND2 side; see `mito-contig-vs-ref.sh`) sits between the two ends of the control region
#' on the circular molecule. A read containing the end of the gene region
#' and then its start (same strand, in molecule order) crossed the whole
#' control region, and the read bases between them give that molecule's
#' control-region length ("wrap" reads). Reads that run past one end of the
#' gene region without wrapping give a lower bound for that molecule.
#' Uses all recruited mitochondrial reads from `mito-reassemble.sh`
#' (`_data/mito-flye-cache/reassemble/mito_reads.fastq`); alignments
#' >= 500 bp (minimap2 map-ont). First checks the wrap detection on
#' synthetic reads (gene-region end + known filler + gene-region start, both
#' strands). Versions are read from `_apptainer/mito-env.yml`. Requires
#' mamba. Run from the repo root. Set ENV_DIR to reuse an environment.
#'
#' Output: `_data/mito-cr-span.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-cr-span.txt
YML=${ROOT}/_bash/_apptainer/mito-env.yml
CACHE=${ROOT}/_data/mito-flye-cache
R1C=${CACHE}/with_opts_assembly.fasta
MR=${CACHE}/reassemble/mito_reads.fastq
THREADS=$(sysctl -n hw.ncpu 2> /dev/null || nproc)
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/env}
cd ${WORK}

[ -s ${R1C} ] && [ -s ${MR} ] \
    || { echo "run _bash/tests/mito-reassemble.sh first" 1>&2; exit 1; }

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

# wraps <PAF> -> read, read length, control-region length
wraps () {
    awk '$11 >= 500 {r = $1; n[r]++; Q[r] = $2; i = n[r]; qs[r, i] = $3
                     qe[r, i] = $4; st[r, i] = $5; ts[r, i] = $8; te[r, i] = $9
                     TL = $7}
         END {for (r in n) {best = -1
                for (a = 1; a <= n[r]; a++) for (b = 1; b <= n[r]; b++) {
                    if (a == b || st[r, a] != st[r, b] || te[r, a] < TL - 1000 \
                        || ts[r, b] > 1000) continue
                    if (st[r, a] == "+" && qe[r, a] <= qs[r, b] + 200)
                        gap = qs[r, b] - qe[r, a] + TL - te[r, a] + ts[r, b]
                    else if (st[r, a] == "-" && qe[r, b] <= qs[r, a] + 200)
                        gap = qs[r, a] - qe[r, b] + TL - te[r, a] + ts[r, b]
                    else continue
                    if (best < 0 || gap < best) best = gap}
                if (best >= 0) print r, Q[r], best}}' $1
}

# synthetic check
python3 -I - genes.fa << 'PYEOF' > syn.fq
import random, sys
random.seed(1)
s = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
rc = lambda x: "".join(comp.get(c, "N") for c in reversed(x.upper()))
for g in (500, 3000, 8000):
    r = s[-3000:] + "".join(random.choice("AATTAATTGC") for _ in range(g)) + s[:3000]
    for lab, seq in (("plus", r), ("minus", rc(r))):
        print("@syn_%s_%d\n%s\n+\n%s" % (lab, g, seq, "5" * len(seq)))
PYEOF
minimap2 -x map-ont genes.fa syn.fq 2> /dev/null > syn.paf
wraps syn.paf | sort > syn.txt

minimap2 -x map-ont -t ${THREADS} genes.fa ${MR} 2> /dev/null \
    | awk '$11 >= 500' > g.paf
wraps g.paf > wrap.txt
# extension beyond each end of the gene region (longest per read and side)
awk '{L = $7; ext = -1}
     $9 >= L - 1000 {ext = ($5 == "+") ? $2 - $4 : $3; side = "past_ND2_end"}
     $8 <= 1000 && ext < 0 {ext = ($5 == "+") ? $3 : $2 - $4; side = "before_rRNA_start"}
     ext >= 0 {print side, $1, ext}' g.paf \
    | sort -k1,1 -k2,2 -k3,3nr | awk '!seen[$1, $2]++' > ext.txt

{
    echo "minimap2 $(minimap2 --version); checked $(date +%Y-%m-%d)"
    echo "gene region: $(region genes) of the round-1 contig (${L} bp)"
    echo "synthetic check (read, detected control-region length; true =" \
        "filler + 0):"
    sed 's/^/  /' syn.txt
    echo "reads: $(awk 'NR % 4 == 2' ${MR} | wc -l | tr -d ' '); with" \
        "alignments >= 500 bp to the gene region:" \
        "$(cut -f 1 g.paf | sort -u | wc -l | tr -d ' ')"
    echo "wrap reads (crossed the whole control region): $(wc -l < wrap.txt | tr -d ' ')"
    [ -s wrap.txt ] && awk '{print $3}' wrap.txt | sort -n \
        | awk '{a[NR] = $1} END {printf "  control-region length: min %d," \
               " median %d, max %d\n", a[1], a[int(NR / 2) + 1], a[NR]}'
    echo "read extension beyond the gene region without wrapping:"
    for s in past_ND2_end before_rRNA_start; do
        awk -v s=${s} '$1 == s {print $3}' ext.txt | sort -n \
            | awk -v s=${s} '{a[NR] = $1; if ($1 >= 5000) c5++
                              if ($1 >= 10000) c10++; if ($1 >= 15000) c15++}
                END {printf "  %-18s reads %5d; median %d, 99th pct %d, max %d;" \
                     " >= 5 kb %d, >= 10 kb %d, >= 15 kb %d\n", s, NR,
                     a[int(NR / 2) + 1], a[int(NR * 0.99) + 1], a[NR], c5, c10, c15}'
    done
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
