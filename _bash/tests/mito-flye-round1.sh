#!/bin/bash

#' Local check (not a cluster job): why does Flye round 1 in
#' `00-references/01-mitogenome.sh` assemble nothing, and does adding
#' `--genome-size 16k --asm-coverage 150` (as in round 2) fix it?
#'
#' Recruits nanopore reads exactly as round 1 does (minimap2 to
#' `tgraci-mito-bait.fasta`, alignments >= 300 bp), then runs Flye with the
#' round-1 command and with the two extra options. Reports, for each run,
#' Flye's exit status, contigs, and whether a 13-22 kb contig exists (the
#' rule `pick_contig` uses). Versions are read from
#' `_apptainer/mito-env.yml`. The read file is streamed, not decompressed.
#' Requires mamba. Takes ~1-2 h. Run from the repo root.
#' Set ENV_DIR to reuse an installed environment; set ONT_FASTQ to change
#' the read file. Recruited reads and assemblies are kept in
#' `_data/mito-flye-cache/` (git-ignored); delete it to re-recruit.
#' `mito-contig-diag.sh` describes the cached assembly.
#'
#' Output: `_data/mito-flye-round1-check.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-flye-round1-check.txt
YML=${ROOT}/_bash/_apptainer/mito-env.yml
BAIT=${ROOT}/_bash/references/tgraci-mito-bait.fasta
COI=${ROOT}/_bash/references/tgraci-coi.fasta
CACHE=${ROOT}/_data/mito-flye-cache
ONT_FASTQ=${ONT_FASTQ:-${HOME}/_data/midgenomes/_sra/basecalls_guppy-5.0.11.fastq.gz}
THREADS=$(sysctl -n hw.ncpu 2> /dev/null || nproc)
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/env}
cd ${WORK}

pin () { grep -E "^\s+- $1=" ${YML} | sed 's/.*=//'; }

if [ ! -d ${ENV_DIR} ]; then
    mamba create -y -q -p ${ENV_DIR} -c conda-forge -c bioconda \
        flye=$(pin flye) minimap2=$(pin minimap2) seqtk=$(pin seqtk) \
        > /dev/null 2>&1 \
        || CONDA_SUBDIR=osx-64 mamba create -y -q -p ${ENV_DIR} \
            -c conda-forge -c bioconda flye=$(pin flye) \
            minimap2=$(pin minimap2) seqtk=$(pin seqtk) > /dev/null 2>&1 \
        || { echo "could not install mito-env tools" 1>&2; exit 1; }
fi
export PATH=${ENV_DIR}/bin:${PATH}

# Round-1 recruitment, as in `recruit` in 01-mitogenome.sh
mkdir -p ${CACHE}
if [ ! -s ${CACHE}/reads1.fastq ]; then
    minimap2 -x map-ont --secondary=no -t ${THREADS} ${BAIT} ${ONT_FASTQ} \
        2> /dev/null \
        | awk -v m=300 '$11 >= m {print $1}' | sort -u > recruited.names
    seqtk subseq ${ONT_FASTQ} recruited.names > ${CACHE}/reads1.fastq
fi
cp ${CACHE}/reads1.fastq ./
N_READS=$(awk 'END {print NR / 4}' reads1.fastq)
N_BP=$(awk 'NR % 4 == 2 {n += length($0)} END {print n}' reads1.fastq)

# summarize <flye dir> <exit status>
summarize () {
    local d=$1 s=$2
    echo "exit status: ${s}"
    grep -h "Overlap-based coverage\|disjointigs\|Total read length\|ERROR" \
        ${d}/flye.log 2> /dev/null | grep -v DEBUG | sed 's/^\[[^]]*\] //'
    if [ -f ${d}/assembly_info.txt ]; then
        echo "contigs (name, length, coverage, circular):"
        awk 'NR > 1 {print "  " $1, $2, $3, $4}' ${d}/assembly_info.txt
        if awk 'NR > 1 && $2 >= 13000 && $2 <= 22000 {f = 1} END {exit !f}' \
            ${d}/assembly_info.txt; then
            echo "13-22 kb contig: yes"
        else
            echo "13-22 kb contig: no"
        fi
        cp ${d}/assembly.fasta ${CACHE}/${d}_assembly.fasta
    else
        echo "no assembly"
    fi
}

flye --nano-raw reads1.fastq --out-dir as_is --threads ${THREADS} \
    > /dev/null 2>&1
S1=$?
flye --nano-raw reads1.fastq --out-dir with_opts --threads ${THREADS} \
    --genome-size 16k --asm-coverage 150 > /dev/null 2>&1
S2=$?

{
    echo "Flye $(flye --version), minimap2 $(minimap2 --version)," \
        "seqtk $(pin seqtk); checked $(date +%Y-%m-%d)"
    echo "reads: ${ONT_FASTQ}"
    echo "recruited (>= 300 bp aligned to tgraci-mito-bait.fasta):" \
        "${N_READS} reads, ${N_BP} bp"
    echo
    echo "== round 1 as in 01-mitogenome.sh"
    summarize as_is ${S1}
    echo
    echo "== round 1 with --genome-size 16k --asm-coverage 150"
    summarize with_opts ${S2}
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
