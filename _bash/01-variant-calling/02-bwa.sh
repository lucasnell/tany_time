#!/bin/bash

#' Align Pool-seq reads to genome assembly using `bwa mem`.
#'
#' It first merges paired-end reads, then separately aligns both merged reads
#' and any reads that can't be merged.
#' It then combines these alignments into a single sorted BAM file
#' that also contains read group info.
#'
#' Reads are mapped competitively to the combined reference from
#' `00-references/02-hologenome` (T. gracilentus plus other species). The output BAM
#' includes reads mapped to any species; `03-mpileup` keeps only
#' T. gracilentus contigs. Reads per species (MAPQ >= 20, primary) are
#' written to `<pool>_species-counts.tsv` (in the output tar.gz) for QC.


. /app/.bashrc
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/l/lnell/dna/trimmed"
export PARENT_DIR_OUT="/staging/l/lnell/dna/bwa"
export HOLO_TAR_FULL_PATH="/staging/l/lnell/hologenome.tar"






#' ========================================================================
#' Inputs
#' ========================================================================

export READ_BASE=$1

require_files ${PARENT_DIR_IN}/trimmed_${READ_BASE}.tar ${HOLO_TAR_FULL_PATH}

export GENOME=hologenome.fasta

export READS1=$(read_tar_name ${PARENT_DIR_IN}/trimmed_${READ_BASE}.tar 1)
check_exit_status "reads file name 1" $?
export READS2=$(read_tar_name ${PARENT_DIR_IN}/trimmed_${READ_BASE}.tar 2)
check_exit_status "reads file name 2" $?
READS1=${READS1%.gz}
READS2=${READS2%.gz}





#' ========================================================================
#' Outputs
#' ========================================================================

# Final files / directories
export OUT_DIR=${READ_BASE}_bwa
export OUT_BAM=${READ_BASE}_bwa.bam
export UNMAPPED_BAM=${READ_BASE}_unmapped.bam
# Intermediates:
export MERGED_READS=merged_trimmed_${READ_BASE}.fastq
export UN_READS1=unmerged_${READS1}
export UN_READS2=unmerged_${READS2}
export MERGED_BAM=${OUT_BAM/.bam/_merged.bam}
export UNMAPPED_MERGED_BAM=${UNMAPPED_BAM/.bam/_merged.bam}
export UN_BAM=${OUT_BAM/.bam/_unmerged.bam}
export UNMAPPED_UN_BAM=${UNMAPPED_BAM/.bam/_unmerged.bam}



# if [ -f ${PARENT_DIR_OUT}/${OUT_DIR}.tar.gz ] && [ -f ${PARENT_DIR_OUT}/${OUT_BAM} ]; then
#     echo "Output files already exist" 1>&2
#     safe_exit 0
# fi


#' ========================================================================
#' Prep for downstream steps.
#' ========================================================================

mkdir ${OUT_DIR}
cd ${OUT_DIR}

tar -xf ${PARENT_DIR_IN}/trimmed_${READ_BASE}.tar -C ./
check_exit_status "extract reads tar file" $?

gunzip ${READS1}.gz && \
    gunzip ${READS2}.gz
check_exit_status "gunzip reads" $?

# Combined reference with its bwa index (built in `00-references/02-hologenome`):
tar -xf ${HOLO_TAR_FULL_PATH} -C ./ \
    && mv hologenome/${GENOME}* ./ \
    && rm -r hologenome
check_exit_status "extract hologenome" $?

# For use later in creating read group info:
export READ_ID=$(head -n 1 ${READS1} | grep -Eo "[ATGCN]+$")


#' ========================================================================
#' Merge paired-end reads.
#' ========================================================================
#'
#' Creating 3 FASTQ files here, for merged reads, unmerged reads # 1,
#' and unmerged reads # 2
N_PAIRS_IN=$(( $(wc -l < ${READS1}) / 4 ))
# (time limit: bbmap 40.02 sometimes hung here with no CPU use)
timeout 6h bbmerge.sh in1=${READS1} in2=${READS2} out=${MERGED_READS} \
    outu1=${UN_READS1} outu2=${UN_READS2} 2> bbmerge.log
check_exit_status "bbmerge" $?
cat bbmerge.log 1>&2
# BBMerge can exit 0 after thread errors that drop most reads (bbmap 40.02,
# 2026-10-29), so check that every pair came out merged or unmerged:
N_PAIRS_OUT=$(( ($(wc -l < ${MERGED_READS}) + $(wc -l < ${UN_READS1})) / 4 ))
echo "bbmerge: ${N_PAIRS_IN} pairs in, ${N_PAIRS_OUT} out (merged + unmerged)"
if grep -q "Exception" bbmerge.log || (( N_PAIRS_IN != N_PAIRS_OUT )); then
    echo "BBMERGE CHECK FAILED: ${N_PAIRS_IN} pairs in, ${N_PAIRS_OUT} out," \
        "$(grep -c Exception bbmerge.log) exceptions" 1>&2
    safe_exit 1
fi
rm bbmerge.log

rm ${READS1} ${READS2}



#' ========================================================================
#' Map unmerged reads.
#' ========================================================================
bwa mem -v 1  -M -t ${THREADS} \
    -R "@RG\tID:${READ_BASE}\tSM:${READ_BASE}_${READ_ID}\tPL:illumina\tLB:lib1" \
    ${GENOME} ${UN_READS1} ${UN_READS2} \
    2> all_${UN_BAM%.bam}.log \
    | samtools view -@ ${THREADS} -bh - \
    > all_${UN_BAM}
check_exit_status "bwa mem (unmerged)" $?
call_bam_stats all_${UN_BAM} "(unmerged, unfiltered)"

rm ${UN_READS1} ${UN_READS2}

samtools view -@ ${THREADS} -bh -f 0x4 all_${UN_BAM} \
    > ${UNMAPPED_UN_BAM}
check_exit_status "samtools view (unmerged, unmapped)" $?
call_bam_stats ${UNMAPPED_UN_BAM} "(unmerged, unmapped)"

samtools view -@ ${THREADS} -bh -q 20 -F 0x100 all_${UN_BAM} \
    > ${UN_BAM}
check_exit_status "samtools view (unmerged, filtered)" $?
call_bam_stats ${UN_BAM} "(unmerged, filtered)"

rm all_${UN_BAM}


#' ========================================================================
#' Map merged reads.
#' ========================================================================
bwa mem -v 1 -M -t ${THREADS} \
    -R "@RG\tID:${READ_BASE}\tSM:${READ_BASE}_${READ_ID}\tPL:illumina\tLB:lib1" \
    ${GENOME} ${MERGED_READS} \
    2> all_${MERGED_BAM%.bam}.log \
    | samtools view -@ ${THREADS} -bh - \
    > all_${MERGED_BAM}
check_exit_status "bwa mem (merged)" $?
call_bam_stats all_${MERGED_BAM} "(merged, unfiltered)"

rm ${MERGED_READS}
rm ${GENOME}*

samtools view -@ ${THREADS} -bh -f 0x4 all_${MERGED_BAM} \
    > ${UNMAPPED_MERGED_BAM}
check_exit_status "samtools view (merged, unmapped)" $?
call_bam_stats ${UNMAPPED_MERGED_BAM} "(merged, unmapped)"

samtools view -@ ${THREADS} -bh -q 20 -F 0x100 all_${MERGED_BAM} \
    > ${MERGED_BAM}
check_exit_status "samtools view (merged, filtered)" $?
call_bam_stats ${MERGED_BAM} "(merged, filtered)"

rm all_${MERGED_BAM}



#' ========================================================================
#' Combine BAM files from mappings of merged and unmerged reads.
#' ========================================================================


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
conda activate gatk-env

gatk MergeSamFiles \
    I=${MERGED_BAM} \
    I=${UN_BAM} \
    SO=coordinate \
    USE_THREADING=true \
    O=${OUT_BAM} \
    VERBOSITY=WARNING
check_exit_status "gatk MergeSamFiles" $?

gatk MergeSamFiles \
    I=${UNMAPPED_MERGED_BAM} \
    I=${UNMAPPED_UN_BAM} \
    SO=coordinate \
    USE_THREADING=true \
    O=${UNMAPPED_BAM} \
    VERBOSITY=WARNING
check_exit_status "Picard_MergeSamFiles (unmapped)" $?

conda deactivate
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


rm ${MERGED_BAM} ${UN_BAM} ${UNMAPPED_MERGED_BAM} ${UNMAPPED_UN_BAM}

call_bam_stats ${OUT_BAM} "(final)"
call_bam_stats ${UNMAPPED_BAM} "(final, unmapped)"


# Reads per species in the final BAM (MAPQ >= 20, primary):
# Groups: "target" = nuclear T. gracilentus; `nt_<code>` = other species'
# nuclear genomes; `mt_Tgraci` = T. gracilentus mitogenome; `mt_<code>` =
# other species' mitogenomes; `bc_<code>` =
# COI barcodes of other taxa (see `00-references/02-hologenome`). Lengths are summed.
samtools index ${OUT_BAM}
check_exit_status "samtools index (final)" $?
samtools idxstats ${OUT_BAM} \
    | awk 'BEGIN{OFS="\t"} $1 != "*" {
        split($1, a, "__"); sp = (a[1] ~ /^(nt|mt|bc)_/) ? a[1] : "target";
        n[sp] += $3; len[sp] += $2 }
        END {print "group", "length", "n_reads";
             for (s in n) print s, len[s], n[s]}' \
    > ${READ_BASE}_species-counts.tsv
check_exit_status "species counts" $?
rm ${OUT_BAM}.bai
cat ${READ_BASE}_species-counts.tsv




#' ========================================================================
#' Handle output files
#' ========================================================================

mv ${OUT_BAM} ${PARENT_DIR_OUT}/

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${PARENT_DIR_OUT}/

rm -r ${OUT_DIR}



safe_exit 0


