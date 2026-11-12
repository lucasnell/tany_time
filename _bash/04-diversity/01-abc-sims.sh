#!/bin/bash

#'
#' ABC test of crashes vs. crashes plus immigration (design:
#' `_R/abc-design.md`): one batch of simulations, or the observed
#' statistics. Both use `bin/abc-stats.R` (statistics) and, for
#' simulations, `bin/abc-sim.py` (msprime, in `sim-env`).
#'
#' Arguments (from `abc_info.txt`, written by `_R/make-job-lists.R`):
#'   1. variant: `main` (reference table), `holdout` (parameter-recovery
#'      checks), or a sensitivity run (`power`, `widen`, `older`, `mu`)
#'   2. model: M1-M5, or `observed`
#'   3. batch number
#'   4. first seed (simulation i uses first seed + i)
#'   5. number of simulations
#'
#' Expected run time per simulation (2,000 loci), projected from 12 prior
#' draws per model with 20 loci on a laptop (2026-11-12): mean 3-10 min by
#' model, up to ~1 h with a large source population (M5 slowest), so ~5-20 h
#' per batch of 100; measured in the pilot (first job per model).
#'


. /app/.bashrc
. bin/helpers.sh || exit 1
conda activate baypass-env

# Edit these for your system:
export INPUTS_TAR="/staging/l/lnell/dna/abc/abc-inputs.tar.gz"
export PARENT_DIR_OUT="/staging/l/lnell/dna/abc/sims"
export ABC_PYTHON=/opt/conda/envs/sim-env/bin/python


export VARIANT=$1
export MODEL=$2
export BATCH=$3
export SEED0=$4
export N_SIMS=$5

if [[ ! "${MODEL}" =~ ^(M[1-5]|observed)$ ]]; then
    echo "Unknown model: ${MODEL}" 1>&2
    safe_exit 1
fi
if [[ ! "${VARIANT}" =~ ^(main|holdout|power|widen|older|mu)$ ]]; then
    echo "Unknown variant: ${VARIANT}" 1>&2
    safe_exit 2
fi

require_files ${INPUTS_TAR}
tar -xzf ${INPUTS_TAR}
check_exit_status "extract inputs" $?

#' Sensitivity runs (`_R/abc-design.md`):
EXTRA=""
case ${VARIANT} in
    power) EXTRA="--ne-form power" ;;
    widen) EXTRA="--widen 10" ;;
    older) EXTRA="--older" ;;
    mu)
        EXTRA="--mu 5.49e-9"
        cp abc-inputs/abc-priors-mu.csv abc-inputs/abc-priors.csv ;;
esac

mkdir -p ${PARENT_DIR_OUT}

if [[ "${MODEL}" == "observed" ]]; then
    OUT=abc-observed-stats.tsv
    Rscript bin/abc-stats.R observed abc-inputs ${OUT}
    check_exit_status "observed stats" $?
    mv ${OUT} ${PARENT_DIR_OUT%/sims}/
    check_exit_status "move output to staging" $?
else
    OUT=abc-${VARIANT}-${MODEL}-b${BATCH}.tsv
    Rscript bin/abc-stats.R sims abc-inputs ${MODEL} ${SEED0} ${N_SIMS} \
        ${OUT} ${EXTRA}
    check_exit_status "simulations" $?
    gzip ${OUT}
    mv ${OUT}.gz ${PARENT_DIR_OUT}/
    check_exit_status "move output to staging" $?
fi

rm -r abc-inputs

safe_exit 0
