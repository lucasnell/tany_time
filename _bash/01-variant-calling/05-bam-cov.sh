#!/bin/bash

#'
#' Calculates coverage in BAM files using `bedtools genomecov`
#'
#' outputs are in BedGraph format, where columns are:
#'
#' chrom chromStart chromEnd dataValue
#'
#' Note that chromosome coordinates are zero-based, half-open.
#'



. /app/.bashrc
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/lnell/dna/bwa"
export PARENT_DIR_OUT="/staging/lnell/dna/bwa"



# export READ_BASE=Ash-19_S5



#' ========================================================================
#' Inputs
#' ========================================================================

export READ_BASE=$1

export IN_BAM=${READ_BASE}_bwa.bam

if [ ! -f ${PARENT_DIR_IN}/${IN_BAM} ]; then
    echo "${PARENT_DIR_IN}/${IN_BAM} does not exist! " 1>&2
    safe_exit 111
fi




#' ========================================================================
#' Outputs
#' ========================================================================


# Final file
export OUT_FILE=${READ_BASE}_bwa_cov.bed.gz




#' ========================================================================
#' Calculate coverage
#' ========================================================================


cp ${PARENT_DIR_IN}/${IN_BAM} ./
check_exit_status "copy BAM" $?

# Takes ~6 min on a 7.5G BAM:
bedtools genomecov -ibam ${IN_BAM} -bga \
    > ${OUT_FILE%.gz}
check_exit_status "calculate coverage" $?

# Takes ~ min on 1.8G BED:
gzip ${OUT_FILE%.gz}
check_exit_status "gzip BED" $?



#' ========================================================================
#' Handle output files
#' ========================================================================

#'
#' Output file will be copied to both staging and (upon exit) ResearchDrive.
#'
#' Copying, not moving so that they will be moved to ResearchDrive after exit
#'
cp ${OUT_FILE} ${PARENT_DIR_OUT}/

rm *.bam *.sif


#' Do NOT use `safe_exit` because that'll remove the files before they
#' can be moved to ResearchDrive.
exit 0
