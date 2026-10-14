#!/bin/bash

#' Basecall one FAST5 chunk (from `01-split-fast5`) on a GPU with Dorado
#' 0.9.6 and the R9.4.1 sup model (`dna_r9.4.1_e8_sup@v3.6`), keeping all
#' reads (quality filtering is in `03-combine`).
#'
#' Runs without the container: the Dorado release bundles its own CUDA
#' libraries and only needs the GPU driver of the execution point, so the
#' helper functions in `_apptainer/helpers.sh` are not available here.
#' Software comes from `00-setup.sh`.
#'
#' Output (`/staging/l/lnell/dna/ont/dorado/`): `<chunk>.bam` (unaligned,
#' with `qs` mean-quality tags) and `<chunk>_summary.tsv`
#' (`dorado summary`).
#'

CHUNK=$1

# Edit these for your system:
export DORADO_DIR="/staging/l/lnell/dorado"
export PARENT_DIR_IN="/staging/l/lnell/dna/ont/fast5_chunks"
export PARENT_DIR_OUT="/staging/l/lnell/dna/ont/dorado"

export DORADO=dorado-0.9.6-linux-x64
export MODEL=dna_r9.4.1_e8_sup@v3.6

fail () {
    echo "Step $1 failed" 1>&2
    cd ${_CONDOR_SCRATCH_DIR:-.} && rm -rf work
    exit ${2:-1}
}

for f in ${DORADO_DIR}/${DORADO}.tar.gz ${DORADO_DIR}/${MODEL}.tar.gz \
    ${PARENT_DIR_IN}/${CHUNK}.tar; do
    [ -f ${f} ] || fail "input check (${f} does not exist)" 111
done

mkdir work && cd work
tar -xzf ${DORADO_DIR}/${DORADO}.tar.gz \
    && tar -xzf ${DORADO_DIR}/${MODEL}.tar.gz \
    || fail "unpack dorado"
mkdir f5 && tar -xf ${PARENT_DIR_IN}/${CHUNK}.tar -C f5 || fail "unpack chunk"
echo "FAST5 files: $(ls f5 | wc -l)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

./${DORADO}/bin/dorado basecaller --device cuda:all ${MODEL} f5/ \
    > ${CHUNK}.bam 2> dorado.log
status=$?
grep -E "Simplex reads basecalled|Samples/s|error|Error" dorado.log
[ ${status} -eq 0 ] || { tail -n 20 dorado.log 1>&2; fail "dorado basecaller"; }

./${DORADO}/bin/dorado summary ${CHUNK}.bam > ${CHUNK}_summary.tsv \
    2>> dorado.log || fail "dorado summary"
echo "reads in summary: $(( $(wc -l < ${CHUNK}_summary.tsv) - 1 ))"

mkdir -p ${PARENT_DIR_OUT}
mv ${CHUNK}.bam ${CHUNK}_summary.tsv ${PARENT_DIR_OUT}/ \
    || fail "move output to staging"

cd ..
rm -rf work
exit 0
