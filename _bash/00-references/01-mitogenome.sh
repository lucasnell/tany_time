#!/bin/bash

#' Assemble the T. gracilentus mitochondrial genome from nanopore reads of a
#' single identified male (R9.4.1 flow cell; re-basecalled with Dorado 0.9.6,
#' model `dna_r9.4.1_e8_sup@v3.6`, in `00-basecall/`).
#'
#' 1. Recruit mitochondrial reads by mapping all reads to T. gracilentus
#'    mitochondrial sequences from NCBI (`tgraci-mito-bait.fasta`).
#' 2. Assemble them with Flye (`--genome-size 16k --asm-coverage 150`;
#'    default settings assemble nothing from these reads,
#'    `_bash/tests/mito-flye-round1.sh`); keep the contig carrying the
#'    T. gracilentus COI barcode MZ659759.1 (most matching bases).
#'    Recruited reads that align better to the nuclear assembly than to
#'    the bait (nuclear copies of mitochondrial DNA, and long nuclear reads
#'    with weak A+T-rich matches) are dropped; see `keep_mito`.
#' 3. Re-recruit reads against that draft (same filter) and re-assemble
#'    with `--genome-size 45k --asm-coverage 150`: these reads cover the
#'    whole molecule, and with `16k` Flye gave only fragments (3 local runs,
#'    and one of two cluster runs), while `45k` gave a 40.3-kb contig that
#'    collapses to 22.8 kb in 3 of 3 local runs
#'    (`_bash/tests/mito-round2-check.sh`).
#' 4. Polish with medaka (model r941_min_sup_g507; medaka has no model for
#'    sup v3.6, and this one did best on the gene region in
#'    `_bash/tests/medaka-model-test.sh`), then with Illumina
#'    reads from one pool (polypolish). Pool reads give a population
#'    consensus, which is what a mapping reference needs.
#' 5. Collapse to one repeat unit (`mito-collapse.sh`): the control region
#'    is > 22 kb and its repeat-unit copy number can't be resolved from
#'    these reads (`_bash/tests/mito-copy-number.sh`; owner's decision
#'    2026-10-26), so the reference is the gene region, the unique
#'    control-region sequence on each side, and one repeat unit.
#' 6. Validate: align published T. gracilentus COI barcodes
#'    (`tgraci-coi.fasta`) to the result. The job fails, without writing
#'    output, unless each barcode in REQUIRED_COI aligns over >= MIN_COI_COV
#'    of its length at >= MIN_COI_ID identity (allowing for COI variation
#'    within species). Only MZ659759.1 is required (owner's decision
#'    2026-10-26): MN679642.1 matches every assembled contig, the pool DNA,
#'    and the adult RNA-seq at only ~94%, so it is reported, not required.
#'    A possible reason (owner, 2026-10-26): MZ659759.1 is from Finland and
#'    MN679642.1 from Canada, so the Canadian barcode may come from a
#'    divergent lineage.
#'
#' Output `Tgraci_mitogenome.tar` (staging) with `Tgraci_mitogenome.fasta`
#' (collapsed), `Tgraci_mitogenome_layout.tsv` (its pieces),
#' `Tgraci_mito_contig.fasta` (the polished contig before collapsing),
#' `coi-check.tsv`, `mt_contigs.txt`, and Flye assembly info.
#'


. /app/.bashrc
. bin/helpers.sh || exit 1
conda activate mito-env

export THREADS=$(count_threads)

# Edit these for your system:
export ONT_FASTQ="/staging/l/lnell/dna/ont/basecalls_dorado-0.9.6_sup.fastq.gz"
export TRIMMED_DIR="/staging/l/lnell/dna/trimmed"
export GENOME_FULL_PATH="/staging/l/lnell/Tgraci_assembly.fasta.gz"
export PARENT_DIR_OUT="/staging/l/lnell"

# Illumina pool used for short-read polishing (high coverage, 80 haploid
# genomes; see `_data/coverage-cutoff-check.csv`):
export POLISH_POOL="SN-1-08_S24"
# Read pairs sampled from that pool for polishing:
export N_PAIRS=4000000
export SEED=2125381094

export MEDAKA_MODEL="r941_min_sup_g507"
export MIN_COI_COV=0.95
export MIN_COI_ID=0.98
export REQUIRED_COI="MZ659759.1"


require_files tgraci-mito-bait.fasta tgraci-coi.fasta \
    tgraci-mito-features.fasta ${ONT_FASTQ} \
    ${TRIMMED_DIR}/trimmed_${POLISH_POOL}.tar ${GENOME_FULL_PATH}

export OUT_DIR=Tgraci_mitogenome
mkdir ${OUT_DIR}
mv tgraci-mito-bait.fasta tgraci-coi.fasta tgraci-mito-features.fasta \
    ${OUT_DIR}/
cd ${OUT_DIR}


#' Choose the assembled mitochondrial contig from Flye output dir $1: the
#' one MZ659759.1 aligns to with the most matching bases (any length; the
#' control-region repeats make the length unreliable). Writes it to file $2.
pick_contig () {
    local contig=$(samtools faidx tgraci-coi.fasta MZ659759.1 \
        | minimap2 -c -x asm20 --secondary=no $1/assembly.fasta - 2> /dev/null \
        | sort -k10,10nr | awk 'NR == 1 {print $6}')
    if [ -z "${contig}" ]; then
        echo "No contig carrying COI MZ659759.1 in $1" 1>&2
        cat $1/assembly_info.txt 1>&2
        safe_exit 1
    fi
    echo "Chose ${contig} from $1:"
    awk -v c=${contig} '$1 == c' $1/assembly_info.txt
    samtools faidx $1/assembly.fasta ${contig} > $2
}

#' Nuclear assembly for `keep_mito`. Its contigs carrying the T. gracilentus
#' COI barcode MZ659759.1 (>= 98% identity over >= 90% of it) are its own
#' mitochondrial contigs (CONTIG34 in the current assembly) and are renamed
#' `mt_<name>`, so they count as mitochondrial; all other contigs count as
#' nuclear.
stage_genome ${GENOME_FULL_PATH}
GENOME=$(basename ${GENOME_FULL_PATH%.gz})
samtools faidx tgraci-coi.fasta MZ659759.1 \
    | minimap2 -c -x asm20 ${GENOME} - 2> /dev/null \
    | awk '$10 / $11 >= 0.98 && $4 - $3 >= 0.9 * $2 {print $6}' \
    | sort -u > mt_contigs.txt
check_exit_status "find mitochondrial contigs in the assembly" $?
echo "Assembly contigs treated as mitochondrial: $(paste -sd ' ' mt_contigs.txt)"
awk 'NR == FNR {mt[$1]; next}
     /^>/ {n = substr($1, 2); $0 = ">" ((n in mt) ? "mt_" : "") n} {print}' \
    mt_contigs.txt ${GENOME} > genome_mt.fasta
rm ${GENOME}

#' Keep reads (names in file $1) whose best alignment (most matching
#' bases; minimap2 map-ont) is to the bait $2 or to a mitochondrial
#' assembly contig rather than to a nuclear contig; overwrites $1.
keep_mito () {
    awk '/^>/ {print ">mt_bait_" substr($1, 2); next} {print}' $2 \
        | cat genome_mt.fasta - > comp.fasta
    seqtk subseq ${ONT_FASTQ} $1 > comp_reads.fastq
    minimap2 -x map-ont --secondary=no -t ${THREADS} comp.fasta \
        comp_reads.fastq 2> /dev/null \
        | awk '!($1 in m) || $10 > m[$1] {m[$1] = $10; t[$1] = $6}
               END {for (r in t) if (t[r] ~ /^mt_/) print r}' \
        | sort > kept.names
    check_exit_status "minimap2 keep_mito" $?
    echo "keep_mito: kept $(wc -l < kept.names) of $(wc -l < $1) recruited reads"
    mv kept.names $1
    rm comp.fasta comp_reads.fastq
}

#' Recruit nanopore reads that align to reference $1 over >= $2 bp, drop
#' those that align better to the nuclear assembly (`keep_mito`); writes
#' them to file $3.
recruit () {
    minimap2 -x map-ont --secondary=no -t ${THREADS} $1 ${ONT_FASTQ} \
        2> /dev/null \
        | awk -v m=$2 '$11 >= m {print $1}' | sort -u > recruited.names
    check_exit_status "minimap2 recruit" $?
    keep_mito recruited.names $1
    seqtk subseq ${ONT_FASTQ} recruited.names > $3
    echo "Recruited $(wc -l < recruited.names) reads to $3"
    rm recruited.names
}



#' ========================================================================
#' Round 1: bait with NCBI sequences
#' ========================================================================

recruit tgraci-mito-bait.fasta 300 reads1.fastq
flye --nano-raw reads1.fastq --out-dir flye1 --threads ${THREADS} \
    --genome-size 16k --asm-coverage 150
check_exit_status "flye round 1" $?
pick_contig flye1 draft1.fasta


#' ========================================================================
#' Round 2: re-recruit with the same-individual draft
#' ========================================================================

recruit draft1.fasta 1000 reads2.fastq
flye --nano-raw reads2.fastq --out-dir flye2 --threads ${THREADS} \
    --genome-size 45k --asm-coverage 150
check_exit_status "flye round 2" $?
cp flye2/assembly_info.txt flye-assembly-info.txt
pick_contig flye2 draft2.fasta


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
    | awk '/^>/ {print ">Tgraci_mito_contig"; next} {print}' \
    > Tgraci_mito_contig.fasta
check_exit_status "polypolish polish" $?


#' ========================================================================
#' Collapse to one repeat unit
#' ========================================================================

mito-collapse.sh Tgraci_mito_contig.fasta tgraci-mito-features.fasta \
    Tgraci_mitogenome
check_exit_status "mito-collapse" $?
echo "Collapsed layout (piece, contig start/end, output start/end, length, truncated):"
cat Tgraci_mitogenome_layout.tsv


#' ========================================================================
#' Validate against published T. gracilentus COI barcodes
#' ========================================================================

#' Columns: barcode, barcode length, aligned barcode bp, identity
#' (matches / alignment length), mitogenome length (collapsed).
MITO_LEN=$(grep -v "^>" Tgraci_mitogenome.fasta | tr -d '\n' | wc -c)
echo -e "barcode\tbarcode_len\taligned_bp\tidentity\tmitogenome_len" \
    > coi-check.tsv
# -c: base-level alignment, so column 10 / 11 is the actual identity
minimap2 -c -x asm20 --secondary=no Tgraci_mitogenome.fasta tgraci-coi.fasta \
    2> /dev/null \
    | awk -v ml=${MITO_LEN} 'BEGIN{OFS="\t"}
        {print $1, $2, $4 - $3, $10 / $11, ml}' >> coi-check.tsv
cat coi-check.tsv

# Required barcodes must align, over enough of their length, at high
# identity (others are only reported in the table):
N_PASS=$(awk -v c=${MIN_COI_COV} -v i=${MIN_COI_ID} -v req="${REQUIRED_COI}" \
    'BEGIN {n = split(req, r, " "); for (k = 1; k <= n; k++) R[r[k]]}
     NR > 1 && ($1 in R) && $3 / $2 >= c && $4 >= i' coi-check.tsv \
    | cut -f 1 | sort -u | wc -l)
N_REQUIRED=$(echo ${REQUIRED_COI} | wc -w)
if (( N_PASS < N_REQUIRED )); then
    echo "COI CHECK FAILED: ${N_PASS} of ${N_REQUIRED} required T. gracilentus" \
        "COI barcodes (${REQUIRED_COI}) align over >= ${MIN_COI_COV} of" \
        "their length at >= ${MIN_COI_ID} identity (see table above)." \
        "No output written." 1>&2
    safe_exit 1
fi
echo "COI CHECK PASSED: ${N_PASS} of ${N_REQUIRED} required barcodes" \
    "(${REQUIRED_COI}); others reported in the table."


#' ========================================================================
#' Outputs
#' ========================================================================

rm -r flye1 flye2 medaka reads*.fastq ill_* medaka.fasta* draft*.fasta* \
    trimmed_* tgraci-mito-bait.fasta tgraci-coi.fasta* genome_mt.fasta \
    tgraci-mito-features.fasta* 2> /dev/null
cd ..
tar -cf ${OUT_DIR}.tar ${OUT_DIR}
mv ${OUT_DIR}.tar ${PARENT_DIR_OUT}/
check_exit_status "move output to staging" $?
rm -r ${OUT_DIR}
if [ -f tany_time.sif ]; then rm tany_time.sif; fi

safe_exit 0
