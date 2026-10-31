#!/bin/bash

#' Local check (not a cluster job): how reliably does Flye round 2 of
#' `00-references/01-mitogenome` give a contig `mito-collapse.sh` can use
#' (full gene region and a full repeat unit)?
#'
#' On 2026-10-27 two cluster runs of the same round 2 gave a 27.1-kb contig
#' (usable) and, the next time, a 6.7-kb contig carrying COI (not usable).
#' This runs Flye as round 2 does (`--nano-raw`, `--genome-size 16k
#' --asm-coverage 150`) REPS times, and with `--genome-size 45k` (as
#' `mito-reassemble.sh`) REPS times, on the filtered mitochondrial reads
#' cached by `mito-reassemble.sh` (`mito_reads.fastq`; close to round 2's
#' reads: whole-molecule bait, nuclear reads dropped). For each run:
#' contigs, which carry COI MZ659759.1, and whether `mito-collapse.sh`
#' succeeds on each contig. Needs ENV_DIR with flye, minimap2, samtools
#' (versions from `_apptainer/mito-env.yml`). Run from the repo root.
#'
#' Output: `_data/mito-round2-check.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-round2-check.txt
MR=${ROOT}/_data/mito-flye-cache/reassemble/mito_reads.fastq
FEATS=${ROOT}/_bash/references/tgraci-mito-features.fasta
COI=${ROOT}/_bash/references/tgraci-coi.fasta
COLLAPSE=${ROOT}/_bash/bin/mito-collapse.sh
REPS=${REPS:-3}
THREADS=$(sysctl -n hw.ncpu 2> /dev/null || nproc)
WORK=$(mktemp -d)
[ -n "${ENV_DIR}" ] && export PATH=${ENV_DIR}/bin:${PATH}
for t in flye minimap2 samtools; do
    command -v ${t} > /dev/null || { echo "${t} not found (set ENV_DIR)" 1>&2; exit 1; }
done
[ -s ${MR} ] || { echo "run _bash/tests/mito-reassemble.sh first" 1>&2; exit 1; }
cd ${WORK}
cp ${FEATS} feats.fa
cp ${COI} coi.fa && samtools faidx coi.fa MZ659759.1 > mz.fa

{
    echo "flye $(flye --version), minimap2 $(minimap2 --version);" \
        "checked $(date +%Y-%m-%d)"
    echo "reads: $(awk 'END {print NR / 4}' ${MR}) (${MR#${ROOT}/})"
    for gs in 16k 45k; do
        for r in $(seq 1 ${REPS}); do
            d=flye_${gs}_${r}
            flye --nano-raw ${MR} --out-dir ${d} --threads ${THREADS} \
                --genome-size ${gs} --asm-coverage 150 > /dev/null 2>&1
            st=$?
            echo
            echo "== --genome-size ${gs}, run ${r} (exit ${st})"
            [ -s ${d}/assembly.fasta ] || { echo "  no assembly"; continue; }
            coi=$(minimap2 -c -x asm20 --secondary=no ${d}/assembly.fasta mz.fa \
                2> /dev/null | sort -k10,10nr | awk 'NR == 1 {print $6}')
            echo "  contig with COI MZ659759.1: ${coi:-none}"
            awk 'NR > 1 {print $1, $2, $3, $4}' ${d}/assembly_info.txt \
                | while read c len cov circ; do
                    samtools faidx ${d}/assembly.fasta ${c} > c.fa
                    if bash ${COLLAPSE} c.fa feats.fa c_mt > /dev/null 2>&1; then
                        res="collapses to $(awk '!/^>/ {n += length($0)} END {print n}' c_mt.fasta) bp"
                    else
                        res="does not collapse"
                    fi
                    printf "  %-10s %6d bp, %6dx, circular %s: %s\n" ${c} ${len} ${cov} ${circ} "${res}"
                done
        done
    done
} > ${OUT}
sed -i '' 's/[[:space:]]*$//' ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
