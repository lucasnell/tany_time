#!/bin/bash

#' Assemble the T. gracilentus mitochondrial genome from nanopore reads of a
#' single identified male (R9.4.1 flow cell; Guppy 5.0.11 hac).
#'
#' 1. Recruit mitochondrial reads by mapping all reads to T. gracilentus
#'    mitochondrial sequences from NCBI (`tgraci-mito-bait.fasta`).
#' 2. Assemble them with Flye; keep the contig 13-22 kb long (preferring
#'    one Flye reports as circular).
#' 3. Re-recruit reads against that draft and re-assemble.
#' 4. Polish with medaka (model r941_min_hac_g507), then with Illumina
#'    reads from one pool (polypolish). Pool reads give a population
#'    consensus, which is what a mapping reference needs.
#' 5. Validate: align published T. gracilentus COI barcodes
#'    (`tgraci-coi.fasta`) to the result. The job fails, without writing
#'    output, unless both align over >= MIN_COI_COV of their length at
#'    >= MIN_COI_ID identity (allowing for COI variation within species).
#'
#' Output `Tgraci_mitogenome.tar` (staging) with
#' `Tgraci_mitogenome.fasta`, `coi-check.tsv`, and Flye assembly info.
#'


. /app/.bashrc
conda activate mito-env

export THREADS=$(count_threads)

# Edit these for your system:
export ONT_FASTQ="/staging/l/lnell/dna/ont/basecalls_guppy-5.0.11.fastq.gz"
export TRIMMED_DIR="/staging/l/lnell/dna/trimmed"
export PARENT_DIR_OUT="/staging/l/lnell"

# Illumina pool used for short-read polishing (high coverage, 80 haploid
# genomes; see `_data/coverage-cutoff-check.csv`):
export POLISH_POOL="SN-1-08_S24"
# Read pairs sampled from that pool for polishing:
export N_PAIRS=4000000
export SEED=2125381094

export MEDAKA_MODEL="r941_min_hac_g507"
export MIN_COI_COV=0.95
export MIN_COI_ID=0.98
export MIN_LEN=13000
export MAX_LEN=22000


require_files tgraci-mito-bait.fasta tgraci-coi.fasta ${ONT_FASTQ} \
    ${TRIMMED_DIR}/trimmed_${POLISH_POOL}.tar

export OUT_DIR=Tgraci_mitogenome
mkdir ${OUT_DIR}
mv tgraci-mito-bait.fasta tgraci-coi.fasta ${OUT_DIR}/
cd ${OUT_DIR}


#' Choose the assembled mitochondrial contig from Flye output dir $1:
#' length in [MIN_LEN, MAX_LEN], circular first, then longest.
#' Writes it to file $2.
pick_contig () {
    local contig=$(awk -v lo=${MIN_LEN} -v hi=${MAX_LEN} \
        'NR > 1 && $2 >= lo && $2 <= hi {print $1, ($4 == "Y" ? 1 : 0), $2}' \
        $1/assembly_info.txt | sort -k2,2nr -k3,3nr | head -n 1 | cut -d' ' -f1)
    if [ -z "${contig}" ]; then
        echo "No contig of ${MIN_LEN}-${MAX_LEN} bp in $1" 1>&2
        cat $1/assembly_info.txt 1>&2
        safe_exit 1
    fi
    echo "Chose ${contig} from $1:"
    awk -v c=${contig} '$1 == c' $1/assembly_info.txt
    samtools faidx $1/assembly.fasta ${contig} > $2
}

#' Recruit nanopore reads that align to reference $1 over >= $2 bp;
#' writes them to file $3.
recruit () {
    minimap2 -x map-ont --secondary=no -t ${THREADS} $1 ${ONT_FASTQ} \
        2> /dev/null \
        | awk -v m=$2 '$11 >= m {print $1}' | sort -u > recruited.names
    check_exit_status "minimap2 recruit" $?
    seqtk subseq ${ONT_FASTQ} recruited.names > $3
    echo "Recruited $(wc -l < recruited.names) reads to $3"
    rm recruited.names
}



#' ========================================================================
#' Round 1: bait with NCBI sequences
#' ========================================================================

recruit tgraci-mito-bait.fasta 300 reads1.fastq
flye --nano-raw reads1.fastq --out-dir flye1 --threads ${THREADS}
check_exit_status "flye round 1" $?
pick_contig flye1 draft1.fasta


#' ========================================================================
#' Round 2: re-recruit with the same-individual draft
#' ========================================================================

recruit draft1.fasta 1000 reads2.fastq
flye --nano-raw reads2.fastq --out-dir flye2 --threads ${THREADS} \
    --genome-size 16k --asm-coverage 150
check_exit_status "flye round 2" $?
pick_contig flye2 draft2.fasta
cp flye2/assembly_info.txt flye-assembly-info.txt


#' ========================================================================
#' Polish: medaka, then polypolish with Illumina reads
#' ========================================================================

medaka_consensus -i reads2.fastq -d draft2.fasta -o medaka \
    -m ${MEDAKA_MODEL} -t ${THREADS}
check_exit_status "medaka" $?

tar -xf ${TRIMMED_DIR}/trimmed_${POLISH_POOL}.tar -C ./
check_exit_status "extract Illumina reads" $?
R1=$(ls trimmed_${POLISH_POOL}*_R1_*.fastq.gz)
R2=$(ls trimmed_${POLISH_POOL}*_R2_*.fastq.gz)
seqtk sample -s ${SEED} ${R1} ${N_PAIRS} > ill_1.fastq \
    && seqtk sample -s ${SEED} ${R2} ${N_PAIRS} > ill_2.fastq
check_exit_status "subsample Illumina reads" $?
rm ${R1} ${R2}

cp medaka/consensus.fasta medaka.fasta
bwa index medaka.fasta
bwa mem -t ${THREADS} -a medaka.fasta ill_1.fastq > ill_1.sam \
    && bwa mem -t ${THREADS} -a medaka.fasta ill_2.fastq > ill_2.sam
check_exit_status "bwa mem (polishing)" $?
polypolish filter --in1 ill_1.sam --in2 ill_2.sam \
    --out1 ill_1.filt.sam --out2 ill_2.filt.sam
check_exit_status "polypolish filter" $?
polypolish polish medaka.fasta ill_1.filt.sam ill_2.filt.sam \
    | awk '/^>/ {print ">Tgraci_mitogenome"; next} {print}' \
    > Tgraci_mitogenome.fasta
check_exit_status "polypolish polish" $?


#' ========================================================================
#' Validate against published T. gracilentus COI barcodes
#' ========================================================================

#' Columns: barcode, barcode length, aligned barcode bp, identity
#' (matches / alignment length), mitogenome length.
MITO_LEN=$(grep -v "^>" Tgraci_mitogenome.fasta | tr -d '\n' | wc -c)
echo -e "barcode\tbarcode_len\taligned_bp\tidentity\tmitogenome_len" \
    > coi-check.tsv
minimap2 -x asm20 --secondary=no Tgraci_mitogenome.fasta tgraci-coi.fasta \
    2> /dev/null \
    | awk -v ml=${MITO_LEN} 'BEGIN{OFS="\t"}
        {print $1, $2, $4 - $3, $10 / $11, ml}' >> coi-check.tsv
cat coi-check.tsv

# Both barcodes must align, over enough of their length, at high identity:
N_PASS=$(awk -v c=${MIN_COI_COV} -v i=${MIN_COI_ID} \
    'NR > 1 && $3 / $2 >= c && $4 >= i' coi-check.tsv | cut -f 1 | sort -u | wc -l)
N_BARCODES=$(grep -c "^>" tgraci-coi.fasta)
if (( N_PASS < N_BARCODES )); then
    echo "COI CHECK FAILED: ${N_PASS} of ${N_BARCODES} T. gracilentus COI" \
        "barcodes align over >= ${MIN_COI_COV} of their length at" \
        ">= ${MIN_COI_ID} identity (see table above). No output written." 1>&2
    safe_exit 1
fi
echo "COI CHECK PASSED: ${N_PASS} of ${N_BARCODES} barcodes."


#' ========================================================================
#' Outputs
#' ========================================================================

rm -r flye1 flye2 medaka reads*.fastq ill_* medaka.fasta* draft*.fasta* \
    trimmed_* tgraci-mito-bait.fasta tgraci-coi.fasta 2> /dev/null
cd ..
tar -cf ${OUT_DIR}.tar ${OUT_DIR}
mv ${OUT_DIR}.tar ${PARENT_DIR_OUT}/
check_exit_status "move output to staging" $?
rm -r ${OUT_DIR}
if [ -f tany_time.sif ]; then rm tany_time.sif; fi

safe_exit 0
