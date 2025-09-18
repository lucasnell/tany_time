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
export INPUTS_TAR_FULL_PATH="/staging/lnell/dna/baypass/baypass-inputs.tar.gz"
export PARENT_DIR_OUT="/staging/lnell/dna/baypass"

#'
#' Note: After checking that the Omega matrices approximately matched
#' between subsamples, I randomly chose to use the third one.
#' This matches the approach from the following papers:
#'   1. Olazcuaga et al., 2020 (doi: 10.1093/molbev/msaa098)
#'   2. Camus et al., 2025 (doi: 10.1101/2024.10.11.617812)
#'
export BP_INIT_TAR_FULL_PATH="/staging/lnell/dna/baypass/baypass-init-sub3.tar.gz"
export OMEGA_FILE="tany-init-sub3_mat_omega.out"


#' ========================================================================
#' Arguments
#' ========================================================================

# export ISING_BETA=1.0
# export RUN_NUM=1
# export SUBSAMPLE=1


export ISING_BETA=$1
export RUN_NUM=$2
export SUBSAMPLE=$3


if (( $(python -c "print(int($ISING_BETA < 0))") )) ||
    (( $(python -c "print(int($ISING_BETA > 1))") )); then
    echo "ISING_BETA outside allowable range: [0,1]" 1>&2
    safe_exit 1
fi
# Only allow two decimals in `ISING_BETA` bc of how `ALL_SEEDS` is generated below
count_decimals () {
    # https://stackoverflow.com/a/61462906/5016095
    decimals=${1#*.}
    echo ${#decimals}
}
if (( $(count_decimals $ISING_BETA) > 2 )); then
    echo "ISING_BETA can have max of 2 decimals." 1>&2
    safe_exit 11
fi

if (( RUN_NUM < 1 )) || (( RUN_NUM > 3 )); then
    echo "RUN_NUM must be >= 1 and <= 3." 1>&2
    echo safe_exit 2
fi
if (( SUBSAMPLE < 1 )) || (( SUBSAMPLE > 5 )); then
    echo "SUBSAMPLE must be >= 1 and <= 5." 1>&2
    echo safe_exit 3
fi

# Seeds for all runs and subsamples (list will vary by `ISING_BETA`):
ALL_SEEDS=($(python3 -c "import random; random.seed(1703756794 + int(100*${ISING_BETA})); print(' '.join([str(random.randint(1, 2147483647)) for x in range(15)]))"))
# Seed for just this run and subsample:
export SEED=${ALL_SEEDS[$(( ($SUBSAMPLE - 1) * 3 + $RUN_NUM - 1 ))]}
if (( SEED < 1 )) || (( SEED > 2147483647 )); then
    echo "SEED too big or small" 1>&2
    safe_exit 4
fi





#' ========================================================================
#' Inputs and output name
#' ========================================================================


if [ ! -f ${INPUTS_TAR_FULL_PATH} ]; then
    echo "${INPUTS_TAR_FULL_PATH} does not exist! " 1>&2
    safe_exit 5
fi
if [ ! -f ${BP_INIT_TAR_FULL_PATH} ]; then
    echo "${BP_INIT_TAR_FULL_PATH} does not exist! " 1>&2
    safe_exit 6
fi


export INPUTS_DIR=$(basename ${INPUTS_TAR_FULL_PATH%.tar.gz})

export OUT_DIR="baypass-beta_${ISING_BETA}-run${RUN_NUM}-sub${SUBSAMPLE}"
export OUT_FILE_PREFIX="tany-beta_${ISING_BETA}-run${RUN_NUM}-sub${SUBSAMPLE}"


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
#' Run baypass (3 independent runs)
#' ========================================================================


export GFILE=$(find ${INPUTS_DIR} -type f -name '*.genobaypass.sub'${SUBSAMPLE})
export SFILE=$(find ${INPUTS_DIR} -type f -name '*.poolsize')
export EFILE=$(find ${INPUTS_DIR} -type f -name '*-log_n.txt')


g_baypass -gfile ${GFILE} \
    -poolsizefile ${SFILE} \
    -efile ${EFILE} \
    -omegafile ${OMEGA_FILE} \
    -auxmodel \
    -isingbeta ${ISING_BETA} \
    -outprefix ${OUT_FILE_PREFIX} \
    -nthreads ${THREADS} \
    -seed ${SEED}
check_exit_status "baypass" $?






#' ========================================================================
#' Deal with outputs
#' ========================================================================


# These aren't necessary to keep for downstream:
rm -r ${INPUTS_DIR} ${OMEGA_FILE}

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

cp ${OUT_DIR}.tar.gz ${PARENT_DIR_OUT}/

rm -r ${OUT_DIR}


exit 0
