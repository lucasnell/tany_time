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
export PARENT_DIR_OUT="/staging/lnell/dna/baypass/beta"

#'
#' Note: Because the Omega matrices based on sub-sampling did not match
#' between subsamples, I am using the one based on all SNPs, despite
#' this not being as computationally efficient.
#'
export BP_INIT_TAR_FULL_PATH="/staging/lnell/dna/baypass/init/baypass-init-sub0.tar.gz"
export OMEGA_FILE="tany-init-sub0_mat_omega.out"


#' ========================================================================
#' Arguments
#' ========================================================================

# export ISING_BETA=1.0
# export RUN_NUM=1


export ISING_BETA=$1
export RUN_NUM=$2


if (( $(python -c "print(int($ISING_BETA < 0))") )) ||
    (( $(python -c "print(int($ISING_BETA > 1))") )); then
    echo "ISING_BETA outside allowable range: [0,1]" 1>&2
    safe_exit 1
fi
if (( RUN_NUM < 1 )) || (( RUN_NUM > 10 )); then
    echo "RUN_NUM must be >= 1 and <= 10." 1>&2
    echo safe_exit 2
fi



if (( SEED < 1 )) || (( SEED > 2147483647 )); then
    echo "SEED must be >= 1 and <= 2147483647" 1>&2
    safe_exit 3
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


export G_FILE=$(find ${INPUTS_DIR} -type f -name '*.genobaypass.sub'${SUBSAMPLE})
export S_FILE=$(find ${INPUTS_DIR} -type f -name '*.poolsize')
export E_FILE=$(find ${INPUTS_DIR} -type f -name '*-log_n.txt')


g_baypass -gfile ${G_FILE} \
    -poolsizefile ${S_FILE} \
    -efile ${E_FILE} \
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
