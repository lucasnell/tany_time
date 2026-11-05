#!/bin/bash

#'
#' Script to use baypass to look for evidence of selection in relation
#' to abundance.
#' Must be run after an initial run of baypass that generates the omega file.
#'


. /app/.bashrc
. bin/helpers.sh || exit 1
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
#' roughly 25–50 h per job if speed per iteration is unchanged. Measured
#' 2026-11-04 for the core model: ~2.5 s per iteration on 1.02 million SNPs
#' (~71 h), so jobs run on thinned subsamples (~176k-341k SNPs; core model
#' measured 2026-11-05: mp050 19-32 h, mp075 18-20 h, mp100 10 h per job,
#' longer than the 13-24 h scaled from SNP counts). The covariate model's speed isn't measured yet, so
#' check the first job's progress before submitting the rest.
export NVAL=5000

#'
#' Each job runs on one thinned subsample (1-3) with that subsample's own
#' Omega from `01-baypass-init` (owner, 2026-11-04); results are combined
#' across subsamples locally. (The 2025 runs used Omega from all SNPs
#' because subsample Omegas disagreed; see `02-check-omegas`.)
#'


#' ========================================================================
#' Arguments
#' ========================================================================

export ISING_BETA=$1
export RUN_NUM=$2
export SEED=$3
#' Missing-data setting tag (e.g., `mp075`; see `00-baypass-inputs.sh`).
export MP_TAG=$4
#' Thinned subsample (1-3); inputs and Ω are those of this setting and
#' subsample.
export SUBSAMPLE=$5
#' Covariate: name of a file in `covariates/` (transferred with the job;
#' written by `_R/baypass-covariates.R`): `hm_logN_6` (primary), `logN`,
#' or `logN_lag1`. One covariate per job: the AUX model with several
#' covariates in one file fits a multiple regression (manual v3.1, 6.2).
export COV=$6
export E_FILE=$(pwd)/${COV}.txt
export INPUTS_TAR_FULL_PATH="${INPUTS_PARENT_DIR}/baypass-inputs_${MP_TAG}_3.tar.gz"
export BP_INIT_TAR_FULL_PATH="${INPUTS_PARENT_DIR}/init/baypass-init-${MP_TAG}-sub${SUBSAMPLE}.tar.gz"
export OMEGA_FILE="tany-init-${MP_TAG}-sub${SUBSAMPLE}_mat_omega.out"


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
if [[ ! "${SUBSAMPLE}" =~ ^[1-3]$ ]]; then
    echo "SUBSAMPLE must be 1, 2, or 3." 1>&2
    safe_exit 4
fi
if [[ ! "${COV}" =~ ^(hm_logN_6|logN|logN_lag1)$ ]]; then
    echo "Unknown covariate: ${COV}" 1>&2
    safe_exit 5
fi





#' ========================================================================
#' Inputs and output name
#' ========================================================================


require_files ${INPUTS_TAR_FULL_PATH} ${BP_INIT_TAR_FULL_PATH} ${E_FILE}


# The `_3` tar holds a folder named without `_3`:
export INPUTS_DIR=$(basename ${INPUTS_TAR_FULL_PATH%_3.tar.gz})

export OUT_DIR="baypass-${MP_TAG}-sub${SUBSAMPLE}-${COV}-beta${ISING_BETA}-run${RUN_NUM}"
export OUT_FILE_PREFIX="tany-${MP_TAG}-sub${SUBSAMPLE}-${COV}-beta${ISING_BETA}-run${RUN_NUM}"


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


export G_FILE=$(find ${INPUTS_DIR} -type f -name '*.genobaypass.sub'${SUBSAMPLE})
export S_FILE=$(find ${INPUTS_DIR} -type f -name '*.poolsize')
# (The inputs' own `tany-log_n.txt` is the same as `covariates/logN.txt`.)
cp ${E_FILE} ./


# Covariates are scaled (mean 0, SD 1) by default since BayPass 2.4
# (`-nocovscaling` turns this off), as the manual recommends (§2.3.3).
g_baypass -gfile ${G_FILE} \
    -poolsizefile ${S_FILE} \
    -efile ${COV}.txt \
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


# Keep this subsample's SNP details (positions) for combining results:
mv ${INPUTS_DIR}/*.snpdet.sub${SUBSAMPLE} ./
check_exit_status "keep SNP details" $?
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

