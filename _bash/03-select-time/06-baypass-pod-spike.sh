#!/bin/bash

#'
#' Power test of candidate statistics (`_R/baypass-pod-spike.R`): the
#' covariate model, as in `03-baypass-beta` and `04-baypass-pod`, on a
#' 100,000-SNP mp075 POD with 1,000 spiked-in associations with
#' `hm_logN_6`, using the mp075 subsample-1 Ω.
#'
#' Arguments (from `baypass_pod_spike_info.txt`, `_R/make-job-lists.R`):
#'   1. run number (1-10)
#'   2. seed
#'
#' Expected run time: ~7-10 h (covariate model, 100,000 SNPs, 8 CPUs;
#' scaled from `04-baypass-pod`).
#'


. /app/.bashrc
. bin/helpers.sh || exit 1
conda activate baypass-env

export THREADS=$(count_threads)

# Edit these for your system:
export INPUTS_PARENT_DIR="/staging/l/lnell/dna/baypass"
export PARENT_DIR_OUT="/staging/l/lnell/dna/baypass/pod"

#' Same as `03-baypass-beta.sh`:
export NVAL=5000

export RUN_NUM=$1
export SEED=$2
export MP_TAG=mp075
export COV=hm_logN_6

export SPIKE_TAR="${INPUTS_PARENT_DIR}/baypass-pod-spike-inputs.tar.gz"
export BP_INIT_TAR_FULL_PATH="${INPUTS_PARENT_DIR}/init/baypass-init-${MP_TAG}-sub1.tar.gz"
export OMEGA_FILE="tany-init-${MP_TAG}-sub1_mat_omega.out"
export POD_DIR="baypass-podspike"
export E_FILE=$(pwd)/${COV}.txt

require_files ${SPIKE_TAR} ${BP_INIT_TAR_FULL_PATH} ${E_FILE}

export OUT_DIR="baypass-podspike-run${RUN_NUM}"
export OUT_FILE_PREFIX="tany-podspike-run${RUN_NUM}"

mkdir ${OUT_DIR}
cd ${OUT_DIR}

tar -xzf ${SPIKE_TAR} ${POD_DIR}
check_exit_status "extract spiked POD" $?
tar -xzf ${BP_INIT_TAR_FULL_PATH} -C ./ \
    && mv ./$(basename ${BP_INIT_TAR_FULL_PATH%.tar.gz})/${OMEGA_FILE} ./ \
    && rm -r $(basename ${BP_INIT_TAR_FULL_PATH%.tar.gz})
check_exit_status "extract omega" $?
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

rm -r ${POD_DIR}
rm -f ${OUT_FILE_PREFIX}_summary_pij.out ${OUT_FILE_PREFIX}_summary_yij.out
gzip ${OUT_FILE_PREFIX}_summary_pi_xtx.out ${OUT_FILE_PREFIX}_summary_betai.out
check_exit_status "gzip" $?

cd ..
tar -cf ${OUT_DIR}.tar ${OUT_DIR}
mkdir -p ${PARENT_DIR_OUT}
mv ${OUT_DIR}.tar ${PARENT_DIR_OUT}/
check_exit_status "move output to staging" $?
rm -r ${OUT_DIR}

safe_exit 0
