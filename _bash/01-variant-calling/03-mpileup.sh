#!/bin/bash

#' This script does the following:
#' - mark and remove duplicates using `Picard MarkDuplicates`
#' - re-align sequences in the proximity of indels with `IndelRealigner` and
#'   `RealignerTargetCreator` in `GATK`
#' - `samtools mpileup`, for T. gracilentus contigs only
#'
#' Steps before mpileup use the combined reference from `00-references/02-hologenome`,
#' matching the alignments from `02-bwa`.



. /app/.bashrc
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/l/lnell/dna/bwa"
export PARENT_DIR_OUT="/staging/l/lnell/dna/mpileup"
export HOLO_TAR_FULL_PATH="/staging/l/lnell/hologenome.tar"





#' ========================================================================
#' Inputs
#' ========================================================================

export READ_BASE=$1

export IN_BAM=${READ_BASE}_bwa.bam
export GENOME=hologenome.fasta

require_files ${PARENT_DIR_IN}/${IN_BAM} ${HOLO_TAR_FULL_PATH}






#' ========================================================================
#' Outputs
#' ========================================================================


# Final files / directories
export OUT_DIR=${READ_BASE}_mpileup
export OUT_FILE=${READ_BASE}_mpileup.txt.gz
# Intermediates:
export MARKDUP_OUT=${IN_BAM/.bam/_nodups.bam}
export REALIGNED_OUT=${MARKDUP_OUT/.bam/_realigned.bam}

# if [ -f ${PARENT_DIR_OUT}/${OUT_DIR}.tar.gz ] && [ -f ${PARENT_DIR_OUT}/${OUT_FILE} ]; then
#     echo "Output files already exist" 1>&2
#     safe_exit 0
# fi


#' ========================================================================
#' Prep for downstream steps.
#' ========================================================================

mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp ${PARENT_DIR_IN}/${IN_BAM} ./
check_exit_status "move BAM" $?

# Combined reference with .fai and .dict (bwa index not needed here):
tar -xf ${HOLO_TAR_FULL_PATH} -C ./ \
    && mv hologenome/${GENOME} hologenome/${GENOME}.fai \
        hologenome/${GENOME%.fasta}.dict hologenome/target_contigs.bed ./ \
    && rm -r hologenome
check_exit_status "extract hologenome" $?



#' ========================================================================
#' Mark and remove duplicates
#' ========================================================================

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
conda activate gatk-env

gatk MarkDuplicates \
    REMOVE_DUPLICATES=true \
    I=${IN_BAM} \
    O=${MARKDUP_OUT} \
    M=${MARKDUP_OUT/.bam/_report.txt} \
    VALIDATION_STRINGENCY=SILENT \
    VERBOSITY=WARNING
check_exit_status "Picard_MarkDuplicates" $?

conda deactivate
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

rm ${IN_BAM}

call_bam_stats ${MARKDUP_OUT} "(nodups)"

samtools index -@ $(($THREADS - 1)) ${MARKDUP_OUT}
check_exit_status "samtools index (nodups)" $?


#' ========================================================================
#' Realign around indels
#' ========================================================================

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
conda activate gatk3-env

GenomeAnalysisTK \
    -T RealignerTargetCreator \
    -nt ${THREADS} \
    -R ${GENOME} \
    -I ${MARKDUP_OUT} \
    -o ${MARKDUP_OUT/.bam/.intervals} \
    -log ${MARKDUP_OUT/.bam/.intervals.log}
check_exit_status "RealignerTargetCreator" $?

# This part can take a while (9+ hours)
GenomeAnalysisTK \
    -T IndelRealigner \
    -R ${GENOME} \
    -I ${MARKDUP_OUT} \
    -targetIntervals ${MARKDUP_OUT/.bam/.intervals} \
    -o ${REALIGNED_OUT}
check_exit_status "IndelRealigner" $?

conda deactivate
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

rm ${MARKDUP_OUT}*

call_bam_stats ${REALIGNED_OUT} "(realigned)"


#' ========================================================================
#' mpileup
#' ========================================================================
# Minimum base quality of 25 matches the threshold `Mpileup2Sync.py` uses
# for allele counts in `02-allele-freqs/01-snape`, so SNAPE frequencies and
# counts come from the same bases (as in DESTv3).
# `-l target_contigs.bed` limits output to T. gracilentus contigs, so reads
# assigned to other species are dropped here.
samtools mpileup -B -Q 25 -l target_contigs.bed -f ${GENOME} \
    ${REALIGNED_OUT} > ${OUT_FILE%.gz}
check_exit_status "samtools mpileup" $?

gzip ${OUT_FILE%.gz}

rm ${GENOME/.fasta/}* target_contigs.bed




#' ========================================================================
#' Handle output files
#' ========================================================================

#'
#' Main mpileup file is separated from directory containing logs.
#' Both the mpileup file and the archived directory will be copied to
#' both staging and (upon exit) ResearchDrive.
#'
mv ${OUT_FILE} ../
cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

# Copying, not moving so that they will be moved to ResearchDrive after exit
cp ${OUT_DIR}.tar.gz ${PARENT_DIR_OUT}/
cp ${OUT_FILE} ${PARENT_DIR_OUT}/

rm -r ${OUT_DIR}

#' Do NOT use `safe_exit` because that'll remove the files before they
#' can be moved to ResearchDrive.
exit 0
