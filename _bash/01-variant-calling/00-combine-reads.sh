#!/bin/bash

#'
#' Convert paired *.fastq.gz files into a single .tar file,
#' then move to staging.
#' The tar step is to reduce the number of files on the cluster.
#'
#' This job can be run with or without the apptainer container.
#'

# Edit this for your system:
export PARENT_DIR_OUT="/staging/lnell/dna/raw_fq"


export READ_BASE=$1
export READS1=${READ_BASE}_L002_R1_001.fastq.gz
export READS2=${READ_BASE}_L002_R2_001.fastq.gz

# Exit with given status but only in non-interactive job:
safe_exit () {
    if [[ $- != *i* ]]; then
        if [ -f tany_time.sif ]; then rm tany_time.sif; fi
        exit $1
    fi
    return 0
}

if [ ! -f ${READS1} ]; then
    echo "${READS1} does not exist! " 1>&2
    safe_exit 1
fi
if [ ! -f ${READS2} ]; then
    echo "${READS2} does not exist! " 1>&2
    safe_exit 1
fi

if [ -f ${PARENT_DIR_OUT}/${READ_BASE}.tar ]; then
    echo "Output file already exists" 1>&2
    safe_exit 0
fi


tar -cf ${READ_BASE}.tar ${READS1} ${READS2}

mv ${READ_BASE}.tar ${PARENT_DIR_OUT}/

rm ${READS1} ${READS2}

safe_exit 0

