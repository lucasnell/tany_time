#!/bin/bash

#' Diagnostic: classify each pool's unmapped reads (reads that mapped to no
#' genome in the combined reference; `02-bwa`) with Kraken2. Not used to
#' include or exclude pools. Two databases:
#'
#' 1. PlusPF-16 (prebuilt, RefSeq-based: archaea, bacteria, viruses,
#'    plasmids, human, UniVec_Core, protozoa, fungi; capped at 16 GB; from
#'    https://benlangmead.github.io/aws-indexes/k2), on a random subsample
#'    of N_READS reads: what the non-target DNA is overall, including
#'    endosymbiotic bacteria (e.g., Wolbachia, Spiroplasma) and viruses.
#' 2. SILVA SSU rRNA (built by `00-references/03-silva-db`), on all unmapped reads:
#'    eukaryotic parasites (fungi, oomycetes, microsporidia, mermithid
#'    nematodes, water mites, ...), whose genomes aren't in RefSeq. rRNA
#'    reads are a small fraction of reads, so no subsampling.
#'
#' One job for all pools so each database is unpacked once.
#'
#' Output: `kraken.tar` with, per pool, `<pool>_kraken-report.txt`
#' (PlusPF) and `<pool>_silva-report.txt` (SILVA).
#'


. /app/.bashrc
conda activate kraken-env

export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/l/lnell/dna/bwa"
export PARENT_DIR_OUT="/staging/l/lnell/dna/kraken"
export KRAKEN_DB_TAR="/staging/l/lnell/kraken2-db/k2_pluspf_16_GB_20260626.tar.gz"
export SILVA_DB_TAR="/staging/l/lnell/kraken2-db/silva-138.2.tar.gz"

# Reads sampled per pool, sampling seed, and Kraken2 confidence threshold
# (fraction of a read's k-mers that must support its assigned taxon):
export N_READS=1000000
export SEED=1490287653
export CONFIDENCE=0.1


require_files all-pools.txt ${KRAKEN_DB_TAR} ${SILVA_DB_TAR}
export POOLS=($(cat all-pools.txt))

export OUT_DIR=kraken
mkdir ${OUT_DIR} db silva_db
tar -xzf ${KRAKEN_DB_TAR} -C db
check_exit_status "extract PlusPF database" $?
tar -xzf ${SILVA_DB_TAR} -C silva_db
check_exit_status "extract SILVA database" $?
# The database files (hash.k2d etc.) may be in a subdirectory:
export DB_DIR=$(dirname $(find db -name hash.k2d | head -n 1))
export SILVA_DIR=$(dirname $(find silva_db -name hash.k2d | head -n 1))
require_files ${DB_DIR}/hash.k2d ${SILVA_DIR}/hash.k2d


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
        > ${pool}_all.fastq
    check_exit_status "unmapped reads to FASTQ (${pool})" $?
    rm -r ${pool}_bwa
    seqtk sample -s ${SEED} ${pool}_all.fastq ${N_READS} > ${pool}_sub.fastq
    check_exit_status "subsample reads (${pool})" $?

    # SILVA, all unmapped reads:
    kraken2 --db ${SILVA_DIR} --threads ${THREADS} \
        --confidence ${CONFIDENCE} \
        --report ${OUT_DIR}/${pool}_silva-report.txt \
        --output /dev/null \
        ${pool}_all.fastq
    check_exit_status "kraken2 SILVA (${pool})" $?
    rm ${pool}_all.fastq

    kraken2 --db ${DB_DIR} --threads ${THREADS} \
        --confidence ${CONFIDENCE} \
        --report ${OUT_DIR}/${pool}_kraken-report.txt \
        --output /dev/null \
        ${pool}_sub.fastq
    check_exit_status "kraken2 (${pool})" $?
    rm ${pool}_sub.fastq

done

rm -r db silva_db all-pools.txt

tar -cf ${OUT_DIR}.tar ${OUT_DIR}
mkdir -p ${PARENT_DIR_OUT}
# Copying, not moving so that it will be moved to ResearchDrive after exit
cp ${OUT_DIR}.tar ${PARENT_DIR_OUT}/
rm -r ${OUT_DIR}
if [ -f tany_time.sif ]; then rm tany_time.sif; fi

#' Do NOT use `safe_exit` because that'll remove the files before they
#' can be moved to ResearchDrive.
exit 0
