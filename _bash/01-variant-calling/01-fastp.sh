#!/bin/bash

#'
#' Use fastp to trim paired-end, Poolseq (DNA) Illumina reads.
#'



. /app/.bashrc
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/lnell/dna/raw_fq"
export PARENT_DIR_OUT="/staging/lnell/dna/trimmed"

export READ_BASE=$1

export READS1=${READ_BASE}_L002_R1_001.fastq.gz
export READS2=${READ_BASE}_L002_R2_001.fastq.gz

status=0
if [ ! -f ${PARENT_DIR_IN}/${READS1} ]; then
    echo "${PARENT_DIR_IN}/${READS1} does not exist! " 1>&2
    status=111
fi
if [ ! -f ${PARENT_DIR_IN}/${READS2} ]; then
    echo "${PARENT_DIR_IN}/${READS2} does not exist! " 1>&2
    status=111
fi
if (( $status != 0 )); then safe_exit $status; fi


# Outputs:
export TRIM_READS1=trimmed_${READS1}
export TRIM_READS2=trimmed_${READS2}
export TRIM_READS_TAR=trimmed_${READ_BASE}.tar
export OUT_DIR=trimmed_${READ_BASE}


mkdir ${OUT_DIR}
cd ${OUT_DIR}

## If you want a progress bar for an interactive job:
# rsync --progress ${PARENT_DIR_IN}/${READS1} ./ &&
#     rsync --progress ${PARENT_DIR_IN}/${READS2} ./

cp ${PARENT_DIR_IN}/${READS1} ./ &&
    cp ${PARENT_DIR_IN}/${READS2} ./
check_exit_status "copy reads files" $?




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



safe_exit 0


