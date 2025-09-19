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

# export SUBSAMPLE=0

export SUBSAMPLE=$1

if (( SUBSAMPLE < 0 )) || (( SUBSAMPLE > 3 )); then
    echo "First arg to 01-baypass-init.sh must be >= 0 and <= 3." 1>&2
    echo safe_exit 1
fi

if (( SUBSAMPLE > 0 )); then
    INPUTS_TAR_FULL_PATH=${INPUTS_TAR_FULL_PATH%.tar.gz}_3.tar.gz
fi

# ALL_SEEDS=($(python3 -c "import random; random.seed(1310488804); print(' '.join([str(random.randint(1, 2147483647)) for x in range(4)]))"))
# Using line below instead of one above in case of differences among systems.
# This is what was output on my computer:
ALL_SEEDS=(115044731 1141272491 1488962297 93102805)
export SEED=${ALL_SEEDS[${SUBSAMPLE}]}

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
if (( SUBSAMPLE > 0 )); then
    INPUTS_DIR=${INPUTS_DIR%_3}
fi

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


if (( SUBSAMPLE > 0 )); then
    G_FILE=$(find ${INPUTS_DIR} -type f -name '*.genobaypass.sub'${SUBSAMPLE})
else
    G_FILE=$(find ${INPUTS_DIR} -type f -name '*.genobaypass')
fi
export G_FILE
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
