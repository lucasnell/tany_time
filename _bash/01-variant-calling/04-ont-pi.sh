#!/bin/bash

#' This is to get a rough estimate of nucleotide diversity and
#' heterozygosity from Nanopore sequencing to use as priors in SNAPE-pooled.


. /app/.bashrc
conda activate main-env


export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/lnell/dna/ont"
export PARENT_DIR_OUT="/staging/lnell/dna/ont"
export GENOME_FULL_PATH="/staging/lnell/Tgraci_assembly.fasta.gz"



# Outputs
export OUT_DIR=ont-pi
export BAM=ont_align.bam
export BCF=ont_variants.bcf
export PI_FILE=ont_variants.windowed.pi

mkdir ${OUT_DIR}
cd ${OUT_DIR}

# Inputs
export ONT_FASTQ=basecalls_guppy-5.0.11.fastq
export GENOME=$(basename ${GENOME_FULL_PATH%.gz})

if [ ! -f ${PARENT_DIR_IN}/${ONT_FASTQ}.gz ]; then
    echo "${PARENT_DIR_IN}/${ONT_FASTQ}.gz does not exist! " 1>&2
    # Don't actually exit if it's an interactive job:
    if [[ $- != *i* ]]; then exit 111; fi
fi
if [ ! -f ${GENOME_FULL_PATH} ]; then
    echo "${GENOME_FULL_PATH} does not exist! " 1>&2
    if [[ $- != *i* ]]; then exit 222; fi
fi



cp ${PARENT_DIR_IN}/${ONT_FASTQ}.gz ./ &&
    gunzip ${ONT_FASTQ}.gz
check_exit_status "move, gunzip ONT reads" $?
cp ${GENOME_FULL_PATH} ./ &&
    gunzip ${GENOME}.gz
check_exit_status "move, gunzip genome" $?

samtools faidx --length 80 ${GENOME}
check_exit_status "index genome" $?


# ------------------------------
# Nanopore reads
# ------------------------------

minimap2 -x map-ont -a -t $((THREADS - 4)) -K 1G -2 ${GENOME} ${ONT_FASTQ} \
    | samtools sort -@ 2 -o ${BAM/.bam/_unfiltered.bam} -
check_exit_status "minimap2" $?

call_bam_stats ${BAM/.bam/_unfiltered.bam} "(unfiltered)"

rm ${ONT_FASTQ}



samtools view -bh -q 20 -F 0x100 --threads 2 ${BAM/.bam/_unfiltered.bam} \
    > ${BAM}
check_exit_status "samtools view (filter)" $?

call_bam_stats ${BAM} "(filtered)"

rm ${BAM/.bam/_unfiltered.bam}


# Takes about an hour
bcftools mpileup -Ou --config ont -f ${GENOME} ${BAM} \
    | bcftools call -P 0.01 -mv --threads 2 -Ob -o ${BCF}
check_exit_status "bcftools mpileup | call" $?

rm ${GENOME}* ${BAM}*


#' Calculate nucleotide diversity in windows
#' This step creates the file ${PI_FILE}

vcftools --bcf ${BCF} --out ${PI_FILE%.windowed.pi} \
    --remove-indels \
    --max-alleles 2 \
    --window-pi 1000 \
    --window-pi-step 1000
check_exit_status "vcftools" $?

rm *.log

#'
#' This returns the mean nucleotide diversity (pi) among all sites.
#' Watterson's estimator of theta is the number of segregating sites divided
#' by the (n-1)th harmonic number (i.e., sum(1/i) for i in 1:(n-1)),
#' where n is the number of haploid samples.
#' Because this is one diploid individual, Watterson's estimator is the
#' same as pi here.
#'
#' The output from this was "0.008763101163069224".
#' I'll use 0.0088 for my prior for SNAPE-pool.
#'
echo $(python3 << EOF
import pandas as pd
win_pi = pd.read_csv("${PI_FILE}", sep="\t")
print(win_pi["PI"].mean())
EOF
)


gzip ${PI_FILE}


# ------------------------------
# Handle output
# ------------------------------

mv ${PI_FILE}.gz ${BCF} ${PARENT_DIR_OUT}/

cd ..
rm -r ${OUT_DIR}



if [ -f tany_time.sif ]; then rm tany_time.sif; fi
