#!/bin/bash

#' Combine partially masked gSYNC files (coverage, repeat, and indel masks
#' only) from all quality-passing pools into one file per sample set
#' (time, space, all) and data type, keeping only sites that SNAPE calls
#' polymorphic (P >= MIN_PP) in at least one pool of the set:
#'
#' - `snape_masked_<set>`: SNAPE-pooled info (allele frequency estimates);
#'   used for allele frequencies in `_R/`.
#' - `count_masked_<set>`: allele counts by nucleotide; used for BayPass
#'   (via poolfstat).
#'
#' Each output tar contains `<name>.sync.gz` (sites passing the SNAPE
#' filter), `<name>.names.gz` (pool names in column order), and
#' `<name>_noblanks.sync.gz` (those sites with data in every pool).
#'
#' Before 2026-10-10, inputs were the fully masked files, so `_noblanks`
#' kept only sites SNAPE called polymorphic in *every* pool. That excluded
#' sites that drifted to fixation in any pool (e.g., in crash years).
#'
#' Pools in each set come from `time-pools.txt` and `space-pools.txt`
#' (transferred from `_bash/samples/`; written by `_R/pool-qc.R` from the
#' QC metrics of `03-pool-qc`).
#'


. /app/.bashrc
conda activate main-env


# Edit this for your system:
export SNAPE_DIR="/staging/lnell/dna/snape"

# Minimum SNAPE P(polymorphic) required in at least one pool:
export MIN_PP=0.9


require_files time-pools.txt space-pools.txt
export TIME_SAMPS=($(cat time-pools.txt))
export SPACE_SAMPS=($(cat space-pools.txt))
export ALL_SAMPS=("${TIME_SAMPS[@]}" "${SPACE_SAMPS[@]}")



#' ========================================================================
#' Inputs
#' ========================================================================

for samp in ${ALL_SAMPS[@]}; do
    require_files ${SNAPE_DIR}/${samp}_snape.tar
done

for samp in ${ALL_SAMPS[@]}; do
    tar -xf ${SNAPE_DIR}/${samp}_snape.tar ${samp}_snape/${samp}_snape_part_masked.sync.gz \
        && mv ${samp}_snape/${samp}_snape_part_masked.sync.gz ./ \
        && rm -r ${samp}_snape \
        && mv ${samp}_snape_part_masked.sync.gz ${samp}.sync.gz
    check_exit_status "cp, rename ${samp}" $?
done



#' ========================================================================
#' Combine files
#' ========================================================================

for kind in snape count; do

    # `combine-syncs.py` keeps SNAPE info with `--snape`, counts without:
    if [ "${kind}" == "snape" ]; then
        COMB_FLAG="--snape"
    else
        COMB_FLAG=""
    fi

    for type in TIME SPACE ALL; do

        #' Using indirect reference (notice `!`) to retrieve one of the
        #' *_SAMPS arrays:
        SAMPS="${type}_SAMPS[@]"
        SAMPS=("${!SAMPS}")

        # Names of output files:
        OUT_BASE=${kind}_masked_$(echo "${type}" | tr '[:upper:]' '[:lower:]')
        OUT=${OUT_BASE}.sync.gz
        OUT_NAMES=${OUT_BASE}.names.gz
        OUT_NOBLANKS=${OUT_BASE}_noblanks.sync.gz

        # Combine sync files (takes ~17 min)
        combine-syncs.py ${COMB_FLAG} --min-pp ${MIN_PP} -o ${OUT%.gz} \
            ${SAMPS[@]/%/.sync.gz}
        check_exit_status "combine-syncs.py ($kind, $type)" $?

        # Create file of names
        for samp in "${SAMPS[@]}"; do
            echo ${samp} >> ${OUT_NAMES%.gz}
        done

        # Remove sites missing data in any pool.
        # This pattern matches both the 6-field (counts) and 7-field
        # (SNAPE) placeholders for missing data.
        cat ${OUT%.gz} | grep -Fv ".:.:.:.:.:." | gzip > ${OUT_NOBLANKS}
        check_exit_status "make_no_blanks ($kind, $type)" $?

        # compress output
        gzip ${OUT_NAMES%.gz} \
            && gzip ${OUT%.gz}
        check_exit_status "gzip out and names ($kind, $type)" $?

        # create output tar file
        # Note that I'm copying (not moving) to staging bc this file will be
        # transferred to ResearchDrive upon exit
        mkdir ${OUT_BASE}
        mv ${OUT} ${OUT_NAMES} ${OUT_NOBLANKS} ./${OUT_BASE}/ \
            && tar -cf ${OUT_BASE}.tar ${OUT_BASE} \
            && cp ${OUT_BASE}.tar ${SNAPE_DIR}/
        check_exit_status "tar, copy over output files ($kind, $type)" $?

        # remove unnecessary folder of outputs
        rm -r ${OUT_BASE}

    done

done



rm *.sync.gz time-pools.txt space-pools.txt
if [ -f tany_time.sif ]; then rm tany_time.sif; fi

#' Do NOT use `safe_exit` because that'll remove the files before they
#' can be moved to ResearchDrive.
exit 0
