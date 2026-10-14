#!/bin/bash

#' Local check (not a cluster job): which change to Flye round 1 in
#' `00-references/01-mitogenome.sh` gives a usable mitogenome after round 2?
#'
#' Round-1 options (all from the cached round-1 reads):
#'   trim:    the `--genome-size 16k --asm-coverage 150` round-1 contig,
#'            cut to its longest stretch with read depth >= 10% of its
#'            maximum (replaces `pick_contig` in round 1)
#'   strict:  only reads aligned >= 1000 bp to the bait (not 300), Flye with
#'            `--genome-size 16k --asm-coverage 150`
#'   meta:    Flye `--meta`
#'   meta_gs: Flye `--meta --genome-size 16k --asm-coverage 150`
#' Untrimmed drafts must pass `pick_contig` (13-22 kb, circular first, then
#' longest). Round 2 is as in the pipeline: reads aligned >= 1000 bp to the
#' draft, Flye `--genome-size 16k --asm-coverage 150`, `pick_contig`.
#' To read the 23.6-GB read file once, round-2 reads are recruited against
#' all drafts together (secondary alignments kept, so a read can count for
#' every draft it aligns to) and then split by draft.
#' Reports each option's drafts and contigs, and COI barcode alignments
#' (minimap2 -c) for the round-2 contig. Polishing is not run, so
#' identities are for unpolished Flye contigs.
#' Requires `_data/mito-flye-cache/` from `mito-flye-round1.sh` and mamba.
#' Takes ~3-4 h. Run from the repo root. Set ENV_DIR to reuse an installed
#' environment; set ONT_FASTQ to change the read file.
#'
#' Output: `_data/mito-round1-options.txt`; drafts and round-2 assemblies
#' in `_data/mito-flye-cache/options/`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-round1-options.txt
YML=${ROOT}/_bash/_apptainer/mito-env.yml
BAIT=${ROOT}/_bash/references/tgraci-mito-bait.fasta
COI=${ROOT}/_bash/references/tgraci-coi.fasta
CACHE=${ROOT}/_data/mito-flye-cache
OPT_DIR=${CACHE}/options
ONT_FASTQ=${ONT_FASTQ:-${HOME}/_data/midgenomes/_sra/basecalls_guppy-5.0.11.fastq.gz}
THREADS=$(sysctl -n hw.ncpu 2> /dev/null || nproc)
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/env}
OPTS=(trim strict meta meta_gs)
cd ${WORK}

if [ ! -s ${CACHE}/reads1.fastq ] || [ ! -s ${CACHE}/with_opts_assembly.fasta ]; then
    echo "run _bash/tests/mito-flye-round1.sh first" 1>&2
    exit 1
fi

pin () { grep -E "^\s+- $1=" ${YML} | sed 's/.*=//'; }

if [ ! -d ${ENV_DIR} ]; then
    PKGS="flye=$(pin flye) minimap2=$(pin minimap2) seqtk=$(pin seqtk)
          samtools=$(pin samtools)"
    mamba create -y -q -p ${ENV_DIR} -c conda-forge -c bioconda ${PKGS} \
        > /dev/null 2>&1 \
        || CONDA_SUBDIR=osx-64 mamba create -y -q -p ${ENV_DIR} \
            -c conda-forge -c bioconda ${PKGS} > /dev/null 2>&1 \
        || { echo "could not install mito-env tools" 1>&2; exit 1; }
fi
export PATH=${ENV_DIR}/bin:${PATH}
mkdir -p ${OPT_DIR}
cp ${CACHE}/reads1.fastq ${CACHE}/with_opts_assembly.fasta ./

# pick_contig as in 01-mitogenome.sh: <flye dir> -> contig name (or empty)
pick () {
    awk 'NR > 1 && $2 >= 13000 && $2 <= 22000 {print $1, ($4 == "Y"), $2}' \
        $1/assembly_info.txt 2> /dev/null \
        | sort -k2,2nr -k3,3nr | head -n 1 | cut -d' ' -f1
}
# write contig $2 of flye dir $1 to $3, renamed to $4
extract () {
    samtools faidx $1/assembly.fasta $2 | sed "1s/.*/>$4/" > $3
}
contigs () {
    awk 'NR > 1 {print "  " $1, $2, $3, $4}' $1/assembly_info.txt 2> /dev/null \
        || echo "  (no assembly)"
}

# per-option notes: note <r1|r2> <option> <text> appends; shown by show
mkdir notes
note () { printf '%s' "$3" >> notes/$1_$2; }
show () { cat notes/$1_$2 2> /dev/null; }

# ---- round 1
# trim
samtools faidx with_opts_assembly.fasta
minimap2 -a -x map-ont --secondary=no -t ${THREADS} with_opts_assembly.fasta \
    reads1.fastq 2> /dev/null | samtools sort -o r1.bam - 2> /dev/null
samtools index r1.bam
REGION=$(samtools depth -a r1.bam \
    | awk 'NR == FNR {if ($3 > m) m = $3; next}
           $3 >= 0.1 * m {if (!run) {s = $2; run = 1}; e = $2;
                          if (e - s > be - bs) {bs = s; be = e; c = $1}; next}
           {run = 0}
           END {print c ":" bs "-" be}' - <(samtools depth -a r1.bam))
samtools faidx with_opts_assembly.fasta ${REGION} | sed "1s/.*/>trim/" \
    > draft_trim.fasta
note r1 trim "core ${REGION} of the 29.5-kb contig"

# strict
minimap2 -x map-ont --secondary=no -t ${THREADS} ${BAIT} reads1.fastq \
    2> /dev/null | awk '$11 >= 1000 {print $1}' | sort -u > strict.names
seqtk subseq reads1.fastq strict.names > reads_strict.fastq
flye --nano-raw reads_strict.fastq --out-dir r1_strict --threads ${THREADS} \
    --genome-size 16k --asm-coverage 150 > /dev/null 2>&1
note r1 strict "$(wc -l < strict.names | tr -d ' ') reads, Flye exit $?"

# meta, meta_gs
flye --nano-raw reads1.fastq --out-dir r1_meta --threads ${THREADS} --meta \
    > /dev/null 2>&1
note r1 meta "Flye exit $?"
flye --nano-raw reads1.fastq --out-dir r1_meta_gs --threads ${THREADS} \
    --meta --genome-size 16k --asm-coverage 150 > /dev/null 2>&1
note r1 meta_gs "Flye exit $?"

for o in strict meta meta_gs; do
    c=$(pick r1_${o})
    if [ -n "${c}" ]; then
        extract r1_${o} ${c} draft_${o}.fasta ${o}
        note r1 ${o} "; pick_contig chose ${c}"
    else
        note r1 ${o} "; pick_contig found no 13-22 kb contig (option stops)"
    fi
done

# ---- round 2: one pass over the read file against all drafts
DRAFTS=$(ls draft_*.fasta 2> /dev/null)
cat ${DRAFTS} > drafts.fasta
cp ${DRAFTS} ${OPT_DIR}/
minimap2 -x map-ont --secondary=yes -N 10 -p 0.5 -t ${THREADS} drafts.fasta \
    ${ONT_FASTQ} 2> /dev/null | awk '$11 >= 1000 {print $6, $1}' \
    | sort -u > r2_hits.txt
cut -d' ' -f2 r2_hits.txt | sort -u > r2_all.names
seqtk subseq ${ONT_FASTQ} r2_all.names > r2_all.fastq

for o in ${OPTS[@]}; do
    [ -s draft_${o}.fasta ] || continue
    awk -v o=${o} '$1 == o {print $2}' r2_hits.txt > r2_${o}.names
    seqtk subseq r2_all.fastq r2_${o}.names > r2_${o}.fastq
    flye --nano-raw r2_${o}.fastq --out-dir r2_${o} --threads ${THREADS} \
        --genome-size 16k --asm-coverage 150 > /dev/null 2>&1
    note r2 ${o} "$(wc -l < r2_${o}.names | tr -d ' ') reads, Flye exit $?"
    cp r2_${o}/assembly.fasta ${OPT_DIR}/r2_${o}_assembly.fasta 2> /dev/null
    c=$(pick r2_${o})
    if [ -n "${c}" ]; then
        extract r2_${o} ${c} r2_${o}_pick.fasta ${o}_r2
        note r2 ${o} "; pick_contig chose ${c}"
    else
        note r2 ${o} "; pick_contig found no 13-22 kb contig"
    fi
done

# ---- report
{
    echo "Flye $(flye --version), minimap2 $(minimap2 --version)," \
        "samtools $(pin samtools), seqtk $(pin seqtk);" \
        "checked $(date +%Y-%m-%d)"
    echo "round-1 reads: $(awk 'END {print NR / 4}' reads1.fastq)" \
        "(cached, >= 300 bp aligned to the bait)"
    for o in ${OPTS[@]}; do
        echo
        echo "== ${o}"
        echo "round 1: $(show r1 ${o})"
        [ "${o}" != "trim" ] && contigs r1_${o}
        [ -s draft_${o}.fasta ] || continue
        echo "round 2: $(show r2 ${o})"
        contigs r2_${o}
        if [ -s r2_${o}_pick.fasta ]; then
            echo "COI on round-2 contig (barcode, barcode bp, start, end," \
                "strand, aligned bp, identity):"
            minimap2 -c -x asm20 --secondary=yes -N 5 r2_${o}_pick.fasta \
                ${COI} 2> /dev/null \
                | awk '{printf "  %s %d %d %d %s %d %.3f\n", $1, $2, $8, $9,
                        $5, $4 - $3, $10 / $11}'
        fi
    done
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
