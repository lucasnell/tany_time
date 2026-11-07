#!/bin/bash

#'
#' BayPass on pseudo-observed data (PODs; manual v3.1, section 4.1) to
#' calibrate the Bayes factor cutoff: PODs are simulated under the core
#' model with no covariate association (`_R/baypass-pod-inputs.R`), so
#' their BFs from the covariate model are a null distribution under the
#' same settings as `03-baypass-beta`.
#'
#' Two kinds of job (4th argument):
#' - a covariate (`hm_logN_6`, `logN`, `logN_lag1`): the covariate model
#'   exactly as in `03-baypass-beta` (`-auxmodel`, Ising beta = 0, the Ω of
#'   subsample 1 of this setting, which the POD was simulated from);
#' - `core`: the core model, to check that Ω re-estimated from the POD
#'   matches the Ω it was simulated from (Gautier 2015), and for XtX.
#'


. /app/.bashrc
. bin/helpers.sh || exit 1
conda activate baypass-env

export THREADS=$(count_threads)

# Edit these for your system:
export INPUTS_PARENT_DIR="/staging/l/lnell/dna/baypass"
export PARENT_DIR_OUT="/staging/l/lnell/dna/baypass/pod"

#' Same as `03-baypass-beta.sh`. Expected run time: ~28 h for 400,000 SNPs
#' (covariate model measured at ~24 h per 341,286 SNPs, 2026-11-06; up to
#' ~46 h on a node as slow as the slowest `01-baypass-init` job).
export NVAL=5000


#' ========================================================================
#' Arguments (from `baypass_pod_info.txt`, written by `_R/make-job-lists.R`)
#' ========================================================================

export RUN_NUM=$1
export SEED=$2
#' Missing-data setting tag (e.g., `mp075`); the POD of this setting.
export MP_TAG=$3
export COV=$4

if (( RUN_NUM < 1 )) || (( RUN_NUM > 100 )); then
    echo "RUN_NUM must be >= 1 and <= 100." 1>&2
    safe_exit 1
fi
if (( SEED < 1 )) || (( SEED > 2147483647 )); then
    echo "SEED must be >= 1 and <= 2147483647" 1>&2
    safe_exit 2
fi
if [[ ! "${MP_TAG}" =~ ^mp(050|075|100)$ ]]; then
    echo "Unknown setting: ${MP_TAG}" 1>&2
    safe_exit 3
fi
if [[ ! "${COV}" =~ ^(hm_logN_6|logN|logN_lag1|core)$ ]]; then
    echo "Unknown covariate: ${COV}" 1>&2
    safe_exit 4
fi


#' ========================================================================
#' Inputs and output name
#' ========================================================================

export POD_TAR_FULL_PATH="${INPUTS_PARENT_DIR}/baypass-pod-inputs.tar.gz"
export BP_INIT_TAR_FULL_PATH="${INPUTS_PARENT_DIR}/init/baypass-init-${MP_TAG}-sub1.tar.gz"
export OMEGA_FILE="tany-init-${MP_TAG}-sub1_mat_omega.out"
export POD_DIR="baypass-pod_${MP_TAG}"

require_files ${POD_TAR_FULL_PATH} ${BP_INIT_TAR_FULL_PATH}
if [[ "${COV}" != "core" ]]; then
    export E_FILE=$(pwd)/${COV}.txt
    require_files ${E_FILE}
fi

export OUT_DIR="baypass-pod-${MP_TAG}-${COV}-run${RUN_NUM}"
export OUT_FILE_PREFIX="tany-pod-${MP_TAG}-${COV}-run${RUN_NUM}"

mkdir ${OUT_DIR}
cd ${OUT_DIR}

tar -xzf ${POD_TAR_FULL_PATH} ${POD_DIR}
check_exit_status "cp, extract POD" $?

# The Ω the POD was simulated from (also given to the core-model job, only
# so it's kept with the output for comparison):
tar -xzf ${BP_INIT_TAR_FULL_PATH} -C ./ \
    && mv ./$(basename ${BP_INIT_TAR_FULL_PATH%.tar.gz})/${OMEGA_FILE} ./ \
    && rm -r $(basename ${BP_INIT_TAR_FULL_PATH%.tar.gz})
check_exit_status "cp, extract omega" $?


#' ========================================================================
#' Run baypass
#' ========================================================================

if [[ "${COV}" == "core" ]]; then
    g_baypass -gfile ${POD_DIR}/tany-pod.genobaypass \
        -poolsizefile ${POD_DIR}/tany-pod.poolsize \
        -nval ${NVAL} \
        -outprefix ${OUT_FILE_PREFIX} \
        -nthreads ${THREADS} \
        -seed ${SEED}
    check_exit_status "baypass core" $?
else
    cp ${E_FILE} ./
    g_baypass -gfile ${POD_DIR}/tany-pod.genobaypass \
        -poolsizefile ${POD_DIR}/tany-pod.poolsize \
        -efile ${COV}.txt \
        -omegafile ${OMEGA_FILE} \
        -auxmodel \
        -isingbeta 0 \
        -nval ${NVAL} \
        -outprefix ${OUT_FILE_PREFIX} \
        -nthreads ${THREADS} \
        -seed ${SEED}
    check_exit_status "baypass" $?
fi


#' ========================================================================
#' Deal with outputs
#' ========================================================================

# Simulated inputs are in `baypass-pod-inputs.tar.gz`; per-pool estimates
# (pij, yij) aren't needed for calibration and are the largest outputs.
rm -r ${POD_DIR}
rm -f ${OUT_FILE_PREFIX}_summary_pij.out ${OUT_FILE_PREFIX}_summary_yij.out
gzip ${OUT_FILE_PREFIX}_summary_pi_xtx.out
check_exit_status "gzip pi_xtx" $?
if [[ "${COV}" != "core" ]]; then
    gzip ${OUT_FILE_PREFIX}_summary_betai.out
    check_exit_status "gzip betai" $?
fi

cd ..
tar -cf ${OUT_DIR}.tar ${OUT_DIR}

mkdir -p ${PARENT_DIR_OUT}
mv ${OUT_DIR}.tar ${PARENT_DIR_OUT}/
check_exit_status "move output to staging" $?

rm -r ${OUT_DIR}


safe_exit 0
