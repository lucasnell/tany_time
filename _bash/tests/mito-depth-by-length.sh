#!/bin/bash

#' Local check (not a cluster job): is nanopore read depth even along the
#' reassembled mitochondrial contig, as a single uniform circular molecule
#' would give?
#'
#' Maps all recruited mitochondrial reads (`mito-reassemble.sh` cache:
#' `_data/mito-flye-cache/reassemble/mito_reads.fastq`) to the new
#' SMARTdenovo contig (`asm_smartdenovo.fa`; gene region once, with unique
#' control-region sequence and then repeat-unit copies on each side) and
#' reports mean depth (primary alignments, MAPQ >= 20) by read-length class
#' in each region. Regions are located on the contig by alignment
#' (minimap2 -c -x asm20) of the round-1 contig's gene region
#' (3600-19500) and one repeat unit (1-2472); the unique control-region
#' segments lie between the gene region and the nearest repeat copy.
#' MAPQ >= 20 undercounts the repeat copies (similar copies give low MAPQ);
#' the unique control-region segments are not affected. Versions are read
#' from `_apptainer/mito-env.yml`. Requires mamba. Run from the repo root.
#' Set ENV_DIR to reuse an environment.
#'
#' Output: `_data/mito-depth-by-length.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-depth-by-length.txt
YML=${ROOT}/_bash/_apptainer/mito-env.yml
CACHE=${ROOT}/_data/mito-flye-cache
R1C=${CACHE}/with_opts_assembly.fasta
ASM=${CACHE}/reassemble/asm_smartdenovo.fa
MR=${CACHE}/reassemble/mito_reads.fastq
THREADS=$(sysctl -n hw.ncpu 2> /dev/null || nproc)
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/env}
cd ${WORK}

[ -s ${ASM} ] && [ -s ${MR} ] \
    || { echo "run _bash/tests/mito-reassemble.sh first" 1>&2; exit 1; }

pin () { grep -E "^\s+- $1=" ${YML} | sed 's/.*=//'; }

if [ ! -d ${ENV_DIR} ]; then
    mamba create -y -q -p ${ENV_DIR} -c conda-forge -c bioconda \
        minimap2=$(pin minimap2) samtools=$(pin samtools) > /dev/null 2>&1 \
        || CONDA_SUBDIR=osx-64 mamba create -y -q -p ${ENV_DIR} \
            -c conda-forge -c bioconda minimap2=$(pin minimap2) \
            samtools=$(pin samtools) > /dev/null 2>&1 \
        || { echo "could not install minimap2 and samtools" 1>&2; exit 1; }
fi
export PATH=${ENV_DIR}/bin:${PATH}

cp ${R1C} r1c.fa
C=$(head -1 r1c.fa | sed 's/^>//; s/ .*//')
samtools faidx r1c.fa ${C}:3600-19500 | sed '1s/.*/>genes/' > genes.fa
samtools faidx r1c.fa ${C}:1-2472 | sed '1s/.*/>unit/' > unit.fa
awk '/^>/ {print ">asm"; next} {print}' ${ASM} > asm.fa

# Regions on the contig: gene region; repeat copies (unit hits >= 1.5 kb);
# unique control-region sequence = between the gene region and the nearest
# repeat copy on each side
minimap2 -c -x asm20 --secondary=no asm.fa genes.fa 2> /dev/null \
    | awk '{print "gene", $8, $9; exit}' > regions.txt
G0=$(awk '{print $2}' regions.txt); G1=$(awk '{print $3}' regions.txt)
minimap2 -c -x asm20 --secondary=yes -N 50 -p 0.3 asm.fa unit.fa 2> /dev/null \
    | awk -v g0=${G0} '$11 >= 1500 {print ($8 < g0 ? "repeats_left" : "repeats_right"), $8, $9}' \
    | sort -k2,2n > units.txt
awk -v g0=${G0} -v g1=${G1} '
    $1 == "repeats_left" && $3 > le {le = $3}
    $1 == "repeats_right" && (!rs || $2 < rs) {rs = $2}
    END {if (le) print "CR_unique_left", le, g0
         if (rs) print "CR_unique_right", g1, rs}' units.txt >> regions.txt
cat units.txt >> regions.txt

minimap2 -a -x map-ont --secondary=no -t ${THREADS} asm.fa ${MR} 2> /dev/null \
    | samtools view -b -F 0x904 -q 20 - | samtools sort -o all.bam - 2> /dev/null
samtools index all.bam

{
    echo "minimap2 $(minimap2 --version), samtools $(pin samtools);" \
        "checked $(date +%Y-%m-%d)"
    echo "contig: $(head -1 ${ASM} | sed 's/^>//'); reads: $(awk 'NR % 4 == 2' ${MR} | wc -l | tr -d ' ')"
    echo "regions (name start end):"
    sed 's/^/  /' regions.txt
    echo "mean depth (primary, MAPQ >= 20) by read length:"
    printf "  %-14s" "read length"
    for r in gene CR_unique_left CR_unique_right repeats_left repeats_right; do
        printf " %16s" ${r}
    done
    printf "\n"
    for cls in "0 5000" "5000 11000" "11000 1000000"; do
        set -- ${cls}
        samtools view -h all.bam | awk -v lo=$1 -v hi=$2 \
            '/^@/ || (length($10) >= lo && length($10) < hi)' \
            | samtools view -b -o c.bam - 2> /dev/null
        samtools index c.bam
        samtools depth -a c.bam > c.depth
        printf "  %-14s" "$1-$2"
        for r in gene CR_unique_left CR_unique_right repeats_left repeats_right; do
            awk -v r=${r} 'NR == FNR {if ($1 == r) {s[++n] = $2; e[n] = $3}; next}
                {for (i = 1; i <= n; i++) if ($2 > s[i] && $2 <= e[i]) {d += $3; k++; break}}
                END {if (k) printf " %16.0f", d / k; else printf " %16s", "NA"}' \
                regions.txt c.depth
        done
        printf "\n"
    done
    echo "nuclear coverage of the whole library: ~250x (mito-long-reads.sh)"
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
