#!/bin/bash

#' This script does the following:
#' - mark and remove duplicates using `Picard MarkDuplicates`
#' - re-align sequences in the proximity of indels with `IndelRealigner` and
#'   `RealignerTargetCreator` in `GATK`
#' - `samtools mpileup`



. /app/.bashrc
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/lnell/dna/bwa"
export PARENT_DIR_OUT="/staging/lnell/dna/mpileup"
export GENOME_FULL_PATH="/staging/lnell/Tgraci_assembly.fasta.gz"



# export READ_BASE=Ash-19_S5
# export READ_BASE=Hrisatjorn_S11



#' ========================================================================
#' Inputs
#' ========================================================================

export READ_BASE=$1

export IN_BAM=${READ_BASE}_bwa.bam
export GENOME=$(basename ${GENOME_FULL_PATH%.gz})

if [ ! -f ${PARENT_DIR_IN}/${IN_BAM} ]; then
    echo "${PARENT_DIR_IN}/${IN_BAM} does not exist! " 1>&2
    # Don't actually exit if it's an interactive job:
    if [[ $- != *i* ]]; then exit 111; fi
fi
if [ ! -f ${GENOME_FULL_PATH} ]; then
    echo "${GENOME_FULL_PATH} does not exist! " 1>&2
    if [[ $- != *i* ]]; then exit 222; fi
fi




#' ========================================================================
#' Outputs
#' ========================================================================


# Final files / directories
export OUT_DIR=${READ_BASE}_mpileup
export OUT_FILE=${READ_BASE}_mpileup.txt.gz
# Intermediates:
export MARKDUP_OUT=${IN_BAM/.bam/_nodups.bam}
export REALIGNED_OUT=${MARKDUP_OUT/.bam/_realigned.bam}

if [ -f ${PARENT_DIR_OUT}/${OUT_DIR}.tar.gz ] && [ -f ${PARENT_DIR_OUT}/${OUT_FILE} ]; then
    echo "Output files already exist" 1>&2
    if [[ $- != *i* ]]; then
        if [ -f tany_time.sif ]; then rm tany_time.sif; fi
        exit 0
    fi
fi


#' ========================================================================
#' Prep for downstream steps.
#' ========================================================================

mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp ${PARENT_DIR_IN}/${IN_BAM} ./
check_exit_status "move BAM" $?

cp ${GENOME_FULL_PATH} ./ && \
    gunzip ${GENOME}.gz
check_exit_status "move, gunzip genome" $?


# Needed downstream
samtools faidx --length 80 ${GENOME}



#' ========================================================================
#' Mark and remove duplicates
#' ========================================================================

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
conda activate gatk-env

gatk CreateSequenceDictionary R=${GENOME}
check_exit_status "Picard_CreateSequenceDictionary" $?

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
samtools mpileup -B -f ${GENOME} ${REALIGNED_OUT} > ${OUT_FILE%.gz}
check_exit_status "samtools mpileup" $?

gzip ${OUT_FILE%.gz}

# Sliding window (500bp with step size of 250bp) of coverage for plotting
window-mpileup.py -r ${GENOME} -s 250 -w 500 ${OUT_FILE}
check_exit_status "window-mpileup.py" $?

# Change default naming to specify window size:
mv ${OUT_FILE/.txt/_window.txt} ${OUT_FILE/.txt/_win500.txt}

rm ${GENOME/.fasta/}*




#' ========================================================================
#' Handle output files
#' ========================================================================


mv ${OUT_FILE} ${PARENT_DIR_OUT}/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${PARENT_DIR_OUT}/

rm -r ${OUT_DIR}


if [ -f tany_time.sif ]; then rm tany_time.sif; fi

