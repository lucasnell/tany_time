#!/bin/bash

#' Convert paired *.fastq.gz files into a single .tar file.
#' I do this to reduce the number of files on the cluster.
#' This file is run interactively.


. /app/.bashrc


for fq1 in *_R1_001.fastq.gz
do
    fq2=${fq1/_R1_/_R2_}
    base=$(echo $fq1 | sed 's/_L002_.*//g')
    tar -cf ${base}.tar ${fq1} ${fq2}
done

