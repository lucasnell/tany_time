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
#' BAMs from `02-bwa` are aligned to the combined reference
#' (`00-references/02-hologenome`); output keeps only T. gracilentus contigs.
#'



. /app/.bashrc
. bin/helpers.sh || exit 1
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/l/lnell/dna/bwa"
export PARENT_DIR_OUT="/staging/l/lnell/dna/bwa"
export HOLO_TAR_FULL_PATH="/staging/l/lnell/hologenome.tar"






#' ========================================================================
#' Inputs
#' ========================================================================

export READ_BASE=$1

export IN_BAM=${READ_BASE}_bwa.bam

require_files ${PARENT_DIR_IN}/${IN_BAM} ${HOLO_TAR_FULL_PATH}




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

tar -xf ${HOLO_TAR_FULL_PATH} -C ./ hologenome/target_contigs.bed \
    && mv hologenome/target_contigs.bed ./ \
    && rm -r hologenome
check_exit_status "extract target contigs" $?

# Takes ~6 min on a 7.5G BAM. Keep only T. gracilentus contigs:
bedtools genomecov -ibam ${IN_BAM} -bga \
    | awk 'NR == FNR {keep[$1] = 1; next} ($1 in keep)' target_contigs.bed - \
    > ${OUT_FILE%.gz}
check_exit_status "calculate coverage" $?
rm target_contigs.bed

# Takes ~ min on 1.8G BED:
gzip ${OUT_FILE%.gz}
check_exit_status "gzip BED" $?



#' ========================================================================
#' Handle output files
#' ========================================================================

mv ${OUT_FILE} ${PARENT_DIR_OUT}/
check_exit_status "move output to staging" $?

rm *.bam *.sif


safe_exit 0
