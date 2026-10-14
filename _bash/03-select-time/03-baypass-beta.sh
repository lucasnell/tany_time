#!/bin/bash

#'
#' Script to use baypass to look for evidence of selection in relation
#' to abundance.
#' Must be run after an initial run of baypass that generates the omega file.
#'


. /app/.bashrc
conda activate baypass-env

export THREADS=$(count_threads)

# Edit these for your system:
export INPUTS_PARENT_DIR="/staging/l/lnell/dna/baypass"
export PARENT_DIR_OUT="/staging/l/lnell/dna/baypass/beta"

#' Number of retained MCMC samples. Set explicitly because BayPass 3.1
#' lowered the default from 1,000 to 250. More samples raise the ceiling on
#' Bayes factors (reached when a SNP's PIP is 1 in every sample) and reduce
#' run-to-run Monte Carlo noise.
#' Run time scales roughly linearly with this.
#' Expected run time: 2.31 runs with 1,000 samples (thin 20, burn-in 5,000,
#' up to 20 pilot runs of 500; <= 35,000 iterations) took 9–12 h with 8
#' CPUs. With 3.1 defaults otherwise (thin 20, burn-in 1,000, 10 pilot runs
#' of 100), 5,000 samples is ~102,000 iterations, i.e., ~3–4x longer:
#' roughly 25–50 h per job if speed per iteration is unchanged (not yet
#' checked for 3.1). CHTC's default run limit is 72 h, so check the first
#' job's progress (log file) before submitting the rest.
export NVAL=5000

#'
#' Note: Because the Omega matrices based on sub-sampling did not match
#' between subsamples, I am using the one based on all SNPs, despite
#' this not being as computationally efficient.
#'


#' ========================================================================
#' Arguments
#' ========================================================================

export ISING_BETA=$1
export RUN_NUM=$2
export SEED=$3
#' Missing-data setting tag (e.g., `mp075`; see `00-baypass-inputs.sh`).
#' Inputs and Ω (from all SNPs, `sub0`) are those of this setting.
export MP_TAG=$4
export INPUTS_TAR_FULL_PATH="${INPUTS_PARENT_DIR}/baypass-inputs_${MP_TAG}.tar.gz"
export BP_INIT_TAR_FULL_PATH="${INPUTS_PARENT_DIR}/init/baypass-init-${MP_TAG}-sub0.tar.gz"
export OMEGA_FILE="tany-init-${MP_TAG}-sub0_mat_omega.out"


if ! awk -v b="${ISING_BETA}" 'BEGIN { exit !(b ~ /^[0-9.]+$/ && b >= 0 && b <= 1) }'; then
    echo "ISING_BETA outside allowable range: [0,1]" 1>&2
    safe_exit 1
fi
if (( RUN_NUM < 1 )) || (( RUN_NUM > 100 )); then
    echo "RUN_NUM must be >= 1 and <= 100." 1>&2
    safe_exit 2
fi
if (( SEED < 1 )) || (( SEED > 2147483647 )); then
    echo "SEED must be >= 1 and <= 2147483647" 1>&2
    safe_exit 3
fi





#' ========================================================================
#' Inputs and output name
#' ========================================================================


require_files ${INPUTS_TAR_FULL_PATH} ${BP_INIT_TAR_FULL_PATH}


export INPUTS_DIR=$(basename ${INPUTS_TAR_FULL_PATH%.tar.gz})

export OUT_DIR="baypass-${MP_TAG}-beta${ISING_BETA}-run${RUN_NUM}"
export OUT_FILE_PREFIX="tany-${MP_TAG}-beta${ISING_BETA}-run${RUN_NUM}"


mkdir ${OUT_DIR}
cd ${OUT_DIR}


# All input files:
tar -xzf ${INPUTS_TAR_FULL_PATH} -C ./
check_exit_status "cp, extract inputs" $?

tar -xzf ${BP_INIT_TAR_FULL_PATH} -C ./ \
    && mv ./$(basename ${BP_INIT_TAR_FULL_PATH%.tar.gz})/${OMEGA_FILE} ./ \
    && rm -r $(basename ${BP_INIT_TAR_FULL_PATH%.tar.gz})
check_exit_status "cp, extract omega" $?




#' ========================================================================
#' Run baypass
#' ========================================================================


export G_FILE=$(find ${INPUTS_DIR} -type f -name '*.genobaypass')
export S_FILE=$(find ${INPUTS_DIR} -type f -name '*.poolsize')
export E_FILE=$(find ${INPUTS_DIR} -type f -name '*-log_n.txt')


# Covariates are scaled (mean 0, SD 1) by default since BayPass 2.4
# (`-nocovscaling` turns this off), as the manual recommends (§2.3.3).
g_baypass -gfile ${G_FILE} \
    -poolsizefile ${S_FILE} \
    -efile ${E_FILE} \
    -omegafile ${OMEGA_FILE} \
    -auxmodel \
    -isingbeta ${ISING_BETA} \
    -nval ${NVAL} \
    -outprefix ${OUT_FILE_PREFIX} \
    -nthreads ${THREADS} \
    -seed ${SEED}
check_exit_status "baypass" $?






#' ========================================================================
#' Deal with outputs
#' ========================================================================


# These aren't necessary to keep for downstream:
rm -r ${INPUTS_DIR} ${OMEGA_FILE}

# These are the largest output files by far (BayPass >= 3.0 splits the
# former `summary_yij_pij.out` into `summary_pij.out` and `summary_yij.out`):
gzip ${OUT_FILE_PREFIX}_summary_pij.out \
    && gzip ${OUT_FILE_PREFIX}_summary_yij.out \
    && gzip ${OUT_FILE_PREFIX}_summary_pi_xtx.out \
    && gzip ${OUT_FILE_PREFIX}_summary_betai.out
check_exit_status "baypass" $?

cd ..
tar -cf ${OUT_DIR}.tar ${OUT_DIR}

mv ${OUT_DIR}.tar ${PARENT_DIR_OUT}/
check_exit_status "move output to staging" $?

rm -r ${OUT_DIR}


safe_exit 0

