#!/bin/bash

#'
#' Test of the bias-corrected POD Ω (`_R/baypass-pod-correct.R`): BayPass
#' core model on 100,000-SNP test PODs simulated from the source Ω (`base`)
#' or the corrected Ω* (`corrected`), to compare the re-estimated Ω with
#' the real-data Ω.
#'
#' Arguments (from `baypass_pod_test_info.txt`, `_R/make-job-lists.R`):
#'   1. setting tag (mp050, mp075, mp100)
#'   2. variant (base, corrected)
#'   3. seed
#'
#' Expected run time: ~5-10 h (core model, 100,000 SNPs, 8 CPUs; scaled
#' from `01-baypass-init`).
#'


. /app/.bashrc
. bin/helpers.sh || exit 1
conda activate baypass-env

export THREADS=$(count_threads)

# Edit these for your system:
export INPUTS_TAR="/staging/l/lnell/dna/baypass/baypass-pod-test-inputs.tar.gz"
export PARENT_DIR_OUT="/staging/l/lnell/dna/baypass/pod"

#' Same as `04-baypass-pod.sh`:
export NVAL=5000

export MP_TAG=$1
export VARIANT=$2
export SEED=$3

if [[ ! "${MP_TAG}" =~ ^mp(050|075|100)$ ]]; then
    echo "Unknown setting: ${MP_TAG}" 1>&2
    safe_exit 1
fi
if [[ ! "${VARIANT}" =~ ^(base|corrected)$ ]]; then
    echo "Unknown variant: ${VARIANT}" 1>&2
    safe_exit 2
fi

require_files ${INPUTS_TAR}

export POD_DIR="baypass-podtest_${MP_TAG}_${VARIANT}"
export OUT_DIR="baypass-podtest-${MP_TAG}-${VARIANT}"
export OUT_FILE_PREFIX="tany-podtest-${MP_TAG}-${VARIANT}"

mkdir ${OUT_DIR}
cd ${OUT_DIR}

tar -xzf ${INPUTS_TAR} ${POD_DIR}
check_exit_status "extract test POD" $?

g_baypass -gfile ${POD_DIR}/tany-pod.genobaypass \
    -poolsizefile ${POD_DIR}/tany-pod.poolsize \
    -nval ${NVAL} \
    -outprefix ${OUT_FILE_PREFIX} \
    -nthreads ${THREADS} \
    -seed ${SEED}
check_exit_status "baypass core" $?

# Keep the Ω the POD was simulated from with the output:
mv ${POD_DIR}/source_mat_omega.out ./
rm -r ${POD_DIR}
rm -f ${OUT_FILE_PREFIX}_summary_pij.out ${OUT_FILE_PREFIX}_summary_yij.out
gzip ${OUT_FILE_PREFIX}_summary_pi_xtx.out

cd ..
tar -cf ${OUT_DIR}.tar ${OUT_DIR}
mkdir -p ${PARENT_DIR_OUT}
mv ${OUT_DIR}.tar ${PARENT_DIR_OUT}/
check_exit_status "move output to staging" $?
rm -r ${OUT_DIR}

safe_exit 0
