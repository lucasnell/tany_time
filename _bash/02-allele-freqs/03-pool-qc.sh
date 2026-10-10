#!/bin/bash

#' Per-pool QC metrics used to decide which pools to include
#' (`_R/pool-qc.R`), computed for all candidate pools in each sample set
#' before any allele-frequency analysis:
#'
#' - From partially masked gSYNC files (`pool-qc-sync.py`), per set:
#'   missing-data rate after masking, number of SNPs, and number of private
#'   SNPs (SNPs found only in that pool among the set's candidates).
#' - From upstream outputs (`pool-qc-files.py`): GC content and duplication
#'   rate from fastp, duplication rate from Picard MarkDuplicates, and
#'   mapping rates from `bamtools stats`.
#'
#' Candidate pools come from `time-pools-candidates.txt` and
#' `space-pools-candidates.txt` (transferred from `_bash/samples/`; written
#' by `_R/make-job-lists.R`).
#'
#' Output: `pool-qc.tar` containing `pool-qc-sync_time.tsv`,
#' `pool-qc-sync_space.tsv`, and `pool-qc-files.tsv`.
#'


. /app/.bashrc
conda activate main-env


# Edit these for your system:
export DNA_DIR="/staging/lnell/dna"
export PARENT_DIR_OUT="/staging/lnell/dna/pool-qc"

# Minimum SNAPE P(polymorphic) for a SNP (same as `01-snape`):
export MIN_PP=0.9


require_files time-pools-candidates.txt space-pools-candidates.txt
export TIME_SAMPS=($(cat time-pools-candidates.txt))
export SPACE_SAMPS=($(cat space-pools-candidates.txt))

export OUT_DIR=pool-qc
mkdir ${OUT_DIR}
mv time-pools-candidates.txt space-pools-candidates.txt ./${OUT_DIR}/
cd ${OUT_DIR}



#' ========================================================================
#' Metrics from upstream outputs (one row per pool)
#' ========================================================================

pool-qc-files.py --header > pool-qc-files.tsv

for samp in "${TIME_SAMPS[@]}" "${SPACE_SAMPS[@]}"; do
    mkdir tmp
    # Missing files are reported as NA by `pool-qc-files.py`, so failures
    # here only print warnings.
    tar -xzf ${DNA_DIR}/trimmed/trimmed_${samp}.tar.gz -C tmp \
        trimmed_${samp}/fastp.json \
        || echo "WARNING: no fastp.json for ${samp}" 1>&2
    tar -xzf ${DNA_DIR}/bwa/${samp}_bwa.tar.gz -C tmp \
        --wildcards "${samp}_bwa/*.stats" \
        || echo "WARNING: no bamtools stats for ${samp}" 1>&2
    tar -xzf ${DNA_DIR}/mpileup/${samp}_mpileup.tar.gz -C tmp \
        ${samp}_mpileup/${samp}_bwa_nodups_report.txt \
        || echo "WARNING: no MarkDuplicates report for ${samp}" 1>&2
    pool-qc-files.py --pool ${samp} \
        --fastp tmp/trimmed_${samp}/fastp.json \
        --bwa-dir tmp/${samp}_bwa \
        --markdup tmp/${samp}_mpileup/${samp}_bwa_nodups_report.txt \
        >> pool-qc-files.tsv
    check_exit_status "pool-qc-files.py (${samp})" $?
    rm -r tmp
done



#' ========================================================================
#' Metrics from partially masked gSYNC files (per sample set)
#' ========================================================================

for type in TIME SPACE; do

    SAMPS="${type}_SAMPS[@]"
    SAMPS=("${!SAMPS}")
    type_lc=$(echo "${type}" | tr '[:upper:]' '[:lower:]')

    # Pools without SNAPE output are skipped (reported in stderr), so they
    # have no row in the output and are excluded by `_R/pool-qc.R`.
    HAVE_SAMPS=()
    for samp in "${SAMPS[@]}"; do
        if [ ! -f ${DNA_DIR}/snape/${samp}_snape.tar ]; then
            echo "WARNING: no SNAPE output for ${samp}; skipping" 1>&2
            continue
        fi
        tar -xf ${DNA_DIR}/snape/${samp}_snape.tar \
            ${samp}_snape/${samp}_snape_part_masked.sync.gz \
            && mv ${samp}_snape/${samp}_snape_part_masked.sync.gz ./${samp}.sync.gz \
            && rm -r ${samp}_snape
        check_exit_status "extract ${samp}" $?
        HAVE_SAMPS+=(${samp})
    done

    pool-qc-sync.py --min-pp ${MIN_PP} -o pool-qc-sync_${type_lc}.tsv \
        ${HAVE_SAMPS[@]/%/.sync.gz}
    check_exit_status "pool-qc-sync.py (${type_lc})" $?

    rm ${HAVE_SAMPS[@]/%/.sync.gz}

done



#' ========================================================================
#' Handle output files
#' ========================================================================

cd ..
tar -cf ${OUT_DIR}.tar ${OUT_DIR}
mkdir -p ${PARENT_DIR_OUT}
# Copying, not moving so that it will be moved to ResearchDrive after exit
cp ${OUT_DIR}.tar ${PARENT_DIR_OUT}/
rm -r ${OUT_DIR}
if [ -f tany_time.sif ]; then rm tany_time.sif; fi

#' Do NOT use `safe_exit` because that'll remove the files before they
#' can be moved to ResearchDrive.
exit 0
