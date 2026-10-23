#!/bin/bash

#' This is to get a rough estimate of nucleotide diversity and
#' heterozygosity from Nanopore sequencing to use as priors in SNAPE-pooled.
#' Reads: the Dorado 0.9.6 sup basecalls (`00-basecall/03-combine`; q >= 10).
#' Outputs are tagged with the basecaller (BASECALLS) so runs on different
#' basecalls don't overwrite each other; the Guppy 5.0.11 hac run wrote
#' `ont_variants.vcf.gz` and `ont_variants.windowed.pi.gz` (untagged).


. /app/.bashrc
conda activate main-env


export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/l/lnell/dna/ont"
export PARENT_DIR_OUT="/staging/l/lnell/dna/ont"
export GENOME_FULL_PATH="/staging/l/lnell/Tgraci_assembly.fasta.gz"



export BASECALLS=dorado-0.9.6_sup

# Outputs
export OUT_DIR=ont-pi
export BAM=ont_align.bam
export VCF=ont_variants_${BASECALLS}.vcf.gz
export PI_FILE=ont_variants_${BASECALLS}.windowed.pi

mkdir ${OUT_DIR}
cd ${OUT_DIR}

# Inputs
export ONT_FASTQ=basecalls_${BASECALLS}.fastq
export GENOME=$(basename ${GENOME_FULL_PATH%.gz})

require_files ${PARENT_DIR_IN}/${ONT_FASTQ}.gz ${GENOME_FULL_PATH}



cp ${PARENT_DIR_IN}/${ONT_FASTQ}.gz ./ &&
    gunzip ${ONT_FASTQ}.gz
check_exit_status "move, gunzip ONT reads" $?
stage_genome ${GENOME_FULL_PATH}

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
# Written as VCF, not BCF: in local tests (2026-10-13), vcftools aborted
# reading several valid toy BCFs (from bcftools 1.22 and 1.24, with vcftools
# 0.1.16 and 0.1.17; trigger not isolated) but read the same calls as VCF.
bcftools mpileup -Ou --config ont -f ${GENOME} ${BAM} \
    | bcftools call -P 0.01 -mv --threads 2 -Oz -o ${VCF}
check_exit_status "bcftools mpileup | call" $?

rm ${GENOME}* ${BAM}*


#' Calculate nucleotide diversity in windows
#' This step creates the file ${PI_FILE}

vcftools --gzvcf ${VCF} --out ${PI_FILE%.windowed.pi} \
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
#' The output from this with the Guppy 5.0.11 hac reads was
#' "0.008763051086537054" (2026-10-23 rerun). I'll use 0.0088 for my prior
#' for SNAPE-pool.
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

mv ${PI_FILE}.gz ${VCF} ${PARENT_DIR_OUT}/

cd ..
rm -r ${OUT_DIR}



safe_exit 0
