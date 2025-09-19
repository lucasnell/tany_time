#!/bin/bash

#'
#' Script to do an initial run to estimate Omega (population covariance matrix).
#'

. /app/.bashrc
conda activate baypass-env

export THREADS=$(count_threads)

# Edit these for your system:
export INPUTS_TAR_FULL_PATH="/staging/lnell/dna/baypass/baypass-inputs.tar.gz"
export PARENT_DIR_OUT="/staging/lnell/dna/baypass/init"


export SUBSAMPLE=$1

if (( SUBSAMPLE < 1 )) || (( SUBSAMPLE > 5 )); then
    echo "First arg to 01-baypass-init.sh must be >= 1 and <= 5." 1>&2
    echo safe_exit 1
fi

# ALL_SEEDS=($(python3 -c "import random; random.seed(1310488804); print(' '.join([str(random.randint(1, 2147483647)) for x in range(5)]))"))
# Using line below instead of one above in case of differences among systems.
# This is what was output on my computer:
ALL_SEEDS=(115044731 1141272491 1488962297 93102805 1099835052)
export SEED=${ALL_SEEDS[$(( $SUBSAMPLE - 1 ))]}

if (( SEED < 1 )) || (( SEED > 2147483647 )); then
    echo "SEED too big or small" 1>&2
    safe_exit 2
fi


#' ========================================================================
#' Inputs and output name
#' ========================================================================


if [ ! -f ${INPUTS_TAR_FULL_PATH} ]; then
    echo "${INPUTS_TAR_FULL_PATH} does not exist! " 1>&2
    safe_exit 3
fi


export INPUTS_DIR=$(basename ${INPUTS_TAR_FULL_PATH%.tar.gz})

export OUT_DIR="baypass-init-sub${SUBSAMPLE}"
export OUT_FILE_PREFIX="tany-init-sub${SUBSAMPLE}"


mkdir ${OUT_DIR}
cd ${OUT_DIR}


# All input files:
tar -xzf ${INPUTS_TAR_FULL_PATH} -C ./
check_exit_status "cp, extract inputs" $?




#' ========================================================================
#' Run baypass
#' ========================================================================



export G_FILE=$(find ${INPUTS_DIR} -type f -name '*.genobaypass.sub'${SUBSAMPLE})
export S_FILE=$(find ${INPUTS_DIR} -type f -name '*.poolsize')

# Initial run to calculate Omega (population covariance matrix)
g_baypass -gfile ${G_FILE} \
    -poolsizefile ${S_FILE} \
    -outprefix ${OUT_FILE_PREFIX} \
    -nthreads ${THREADS} \
    -seed ${SEED}
check_exit_status "baypass init" $?




#' ========================================================================
#' Deal with outputs
#' ========================================================================

# These aren't necessary to keep for downstream:
rm -r ${INPUTS_DIR}

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

cp ${OUT_DIR}.tar.gz ${PARENT_DIR_OUT}/

rm -r ${OUT_DIR}


exit 0
