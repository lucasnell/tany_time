#!/bin/bash

#' Local check (not a cluster job): does `locate_features`
#' (`mito-features.sh`) recover the fixed coordinates the mito checks used
#' on the Guppy round-1 contig, and where are the features on the current
#' (sup) round-1 contig?
#'
#' Needs minimap2 and samtools: set ENV_DIR to an environment with both
#' (e.g., the one `medaka-model-test.sh` builds), or they must be on PATH.
#' Run from the repo root.
#'
#' Output: `_data/mito-features-check.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-features-check.txt
GUPPY=${ROOT}/_data/mito-flye-cache-guppy/with_opts_assembly.fasta
SUP=${ROOT}/_data/mito-flye-cache/with_opts_assembly.fasta
EXPECTED="CR_left 1 3600
CR_right 19500 29542
ND2_side 17400 19500
core 9674 15798
genes 3600 19500
rRNA_side 3700 6500
unit 1 2472"
export LC_ALL=C  # sort order of EXPECTED
WORK=$(mktemp -d)
[ -n "${ENV_DIR}" ] && export PATH=${ENV_DIR}/bin:${PATH}
for t in minimap2 samtools; do
    command -v ${t} > /dev/null || { echo "${t} not found (set ENV_DIR)" 1>&2; exit 1; }
done
cd ${WORK}

. ${ROOT}/_bash/tests/mito-features.sh

{
    echo "minimap2 $(minimap2 --version); checked $(date +%Y-%m-%d)"
    for x in "Guppy round-1 contig:${GUPPY}" "sup round-1 contig:${SUP}"; do
        f=${x#*:}
        echo
        echo "== ${x%%:*} (${f#${ROOT}/})"
        locate_features ${f} 2>&1 || continue
        echo "length ${L}; orientation $(cmp -s <(grep -v '^>' r1c.fa | tr -d '\n') \
            <(awk '/^>/ && n++ {exit} !/^>/' ${f} | tr -d '\n') \
            && echo as-is || echo reverse-complemented)"
        sort features.txt | sed 's/^/  /'
        if [ "${f}" == "${GUPPY}" ]; then
            if [ "$(sort features.txt)" == "${EXPECTED}" ]; then
                echo "matches the fixed coordinates used before 2026-10-23: PASS"
            else
                echo "differs from the fixed coordinates used before 2026-10-23: FAIL"
            fi
        fi
    done
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
