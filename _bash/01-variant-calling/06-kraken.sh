#!/bin/bash

#' Diagnostic: classify a random subsample of each pool's unmapped reads
#' (reads that mapped to no genome in the combined reference; `02-bwa`)
#' with Kraken2, to see what the non-target DNA is (e.g., bacteria, fungi,
#' human). Not used to include or exclude pools.
#'
#' One job for all pools so the database is unpacked once.
#'
#' Database: a prebuilt Kraken2 index from
#' https://benlangmead.github.io/aws-indexes/k2 (PlusPF-16: archaea,
#' bacteria, viruses, plasmids, human, UniVec_Core, protozoa, fungi;
#' capped at 16 GB). Download the .tar.gz to `KRAKEN_DB_TAR` on staging.
#'
#' Output: `kraken.tar` with one Kraken2 report per pool
#' (`<pool>_kraken-report.txt`).
#'


. /app/.bashrc
conda activate kraken-env

export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/lnell/dna/bwa"
export PARENT_DIR_OUT="/staging/lnell/dna/kraken"
export KRAKEN_DB_TAR="/staging/lnell/kraken2-db/k2_pluspf_16gb.tar.gz"

# Reads sampled per pool, sampling seed, and Kraken2 confidence threshold
# (fraction of a read's k-mers that must support its assigned taxon):
export N_READS=1000000
export SEED=1490287653
export CONFIDENCE=0.1


require_files all-pools.txt ${KRAKEN_DB_TAR}
export POOLS=($(cat all-pools.txt))

export OUT_DIR=kraken
mkdir ${OUT_DIR} db
tar -xzf ${KRAKEN_DB_TAR} -C db
check_exit_status "extract Kraken2 database" $?
# The database files (hash.k2d etc.) may be in a subdirectory:
export DB_DIR=$(dirname $(find db -name hash.k2d | head -n 1))
require_files ${DB_DIR}/hash.k2d


for pool in ${POOLS[@]}; do

    if [ ! -f ${PARENT_DIR_IN}/${pool}_bwa.tar.gz ]; then
        echo "WARNING: no bwa output for ${pool}; skipping" 1>&2
        continue
    fi

    tar -xzf ${PARENT_DIR_IN}/${pool}_bwa.tar.gz \
        ${pool}_bwa/${pool}_unmapped.bam
    check_exit_status "extract unmapped reads (${pool})" $?

    # Reads are treated as single reads (merged and unmerged mates alike).
    samtools fastq -@ 2 ${pool}_bwa/${pool}_unmapped.bam 2> /dev/null \
        | seqtk sample -s ${SEED} - ${N_READS} \
        > ${pool}_sub.fastq
    check_exit_status "subsample reads (${pool})" $?
    rm -r ${pool}_bwa

    kraken2 --db ${DB_DIR} --threads ${THREADS} \
        --confidence ${CONFIDENCE} \
        --report ${OUT_DIR}/${pool}_kraken-report.txt \
        --output /dev/null \
        ${pool}_sub.fastq
    check_exit_status "kraken2 (${pool})" $?
    rm ${pool}_sub.fastq

done

rm -r db all-pools.txt

tar -cf ${OUT_DIR}.tar ${OUT_DIR}
mkdir -p ${PARENT_DIR_OUT}
# Copying, not moving so that it will be moved to ResearchDrive after exit
cp ${OUT_DIR}.tar ${PARENT_DIR_OUT}/
rm -r ${OUT_DIR}
if [ -f tany_time.sif ]; then rm tany_time.sif; fi

#' Do NOT use `safe_exit` because that'll remove the files before they
#' can be moved to ResearchDrive.
exit 0
