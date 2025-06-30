#!/bin/bash

# Use fastp to trim paired-end, Poolseq (DNA) Illumina reads.



. /app/.bashrc
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/lnell/dna/raw_fq"
export PARENT_DIR_OUT="/staging/lnell/dna/trimmed"



export READ_BASE=$1


if [ ! -f ${PARENT_DIR_IN}/${READ_BASE}.tar ]; then
    echo "${PARENT_DIR_IN}/${READ_BASE}.tar does not exist! " 1>&2
    # Don't actually exit if it's an interactive job:
    if [[ $- != *i* ]]; then
        exit 111
    fi
fi


export READS1=$(read_tar_name ${PARENT_DIR_IN}/${READ_BASE}.tar 1)
check_exit_status "reads file name 1" $?
export READS2=$(read_tar_name ${PARENT_DIR_IN}/${READ_BASE}.tar 2)
check_exit_status "reads file name 2" $?


# Outputs:
export TRIM_READS1=trimmed_${READS1}
export TRIM_READS2=trimmed_${READS2}
export TRIM_READS_TAR=trimmed_${READ_BASE}.tar
export OUT_DIR=trimmed_${READ_BASE}

if [ -f ${PARENT_DIR_OUT}/${OUT_DIR}.tar.gz ] && [ -f ${PARENT_DIR_OUT}/${TRIM_READS_TAR} ]; then
    echo "Output files already exist" 1>&2
    if [[ $- != *i* ]]; then
        if [ -f tany_time.sif ]; then rm tany_time.sif; fi
        exit 0
    fi
fi



mkdir ${OUT_DIR}
cd ${OUT_DIR}

## If you want a progress bar for an interactive job:
# pv ${PARENT_DIR_IN}/${READ_BASE}.tar | tar -x -C ./

tar -xf ${PARENT_DIR_IN}/${READ_BASE}.tar -C ./
check_exit_status "extract reads tar file" $?




# The main things happening here are...
# Automatic adapter trimming (on by default)
# polyG tail trimming for NovaSeq sequencing (on by default)
# enable base correction in overlapped regions for PE data (`--correction`)
# no quality filtering (`--disable_quality_filtering`)
#
# I'm not quality-trimming because `bwa-mem` will soft-mask low-quality reads
# during the alignment phase.

fastp --in1 ${READS1} --in2 ${READS2} \
    --out1 ${TRIM_READS1} --out2 ${TRIM_READS2} \
    --thread ${THREADS} \
    --correction \
    --disable_quality_filtering
check_exit_status "fastp" $?

rm ${READS1} ${READS2}

mkdir ${READ_BASE}_fastqc

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
conda activate qc-env
fastqc ${TRIM_READS1} ${TRIM_READS2} -o ${READ_BASE}_fastqc
check_exit_status "fastqc" $?
conda deactivate
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


tar -cf ${TRIM_READS_TAR} ${TRIM_READS1} ${TRIM_READS2} && \
    mv ${TRIM_READS_TAR} ${PARENT_DIR_OUT}/
check_exit_status "tar, move output reads" $?

rm ${TRIM_READS1} ${TRIM_READS2}


cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR} && \
    mv ${OUT_DIR}.tar.gz ${PARENT_DIR_OUT}/
check_exit_status "tar, move output folder" $?

rm -r ${OUT_DIR}



if [ -f tany_time.sif ]; then rm tany_time.sif; fi
