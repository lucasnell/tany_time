#!/bin/bash

#' Local check (not a cluster job): reassemble the T. gracilentus
#' mitogenome from all mitochondrial nanopore reads with four assemblers,
#' and test each assembly's control-region repeat structure.
#'
#' 1. Recruit reads aligned >= 1 kb (minimap2 map-ont) to a bait of all
#'    mitochondrial contigs seen so far: the cached round-1 contig and the
#'    2019 NextDenovo (ctg000860) and SMARTdenovo (utg319) contigs, which
#'    together carry the gene region, the control-region ends, and copies
#'    of the ~2.5-kb repeat unit (`mito-old-assemblies.sh`). Unlike round 1
#'    of `01-mitogenome` (COI bait only), this samples the whole molecule.
#' 2. Keep the longest, best reads up to TARGET_BP (Filtlong).
#' 3. Assemble with Flye (`--genome-size 45k --asm-coverage 200`), Raven,
#'    NextDenovo (genome_size 45k), and SMARTdenovo.
#' 4. For each contig with >= 5 kb of the gene region (plus the 2019
#'    NextDenovo contig for comparison), report the layout: gene-region
#'    copies, control-region pieces, repeat-unit copies (hits of one unit,
#'    from the round-1 contig), and the distance between gene-region starts.
#'    Pieces of the round-1 contig (gene region, control-region sides, one
#'    unit) are located by alignment (`mito-features.sh`; on the 29.5-kb
#'    Guppy contig 3600-19500, 1-3600, 19500-end, and 1-2472).
#' 5. Test the repeat structure:
#'    a. bridging: recruited reads (MAPQ >= 20, primary) whose single
#'       alignment spans 1 kb either side of each boundary between
#'       neighbouring repeat units; a boundary no read spans is unsupported
#'    b. copy number: Illumina DNA (first ILL_PAIRS pairs of pool
#'       MyKS-19-B_S18) mean depth (all mapping qualities, so reads that
#'       fit several copies still count once) on each unit relative to the
#'       gene region, and their sum: the number of unit copies in the
#'       genome's molecules if Illumina covers the A+T-rich units as well as
#'       the gene region (it tends to under-cover A+T-rich sequence, so the
#'       sum is more likely low than high)
#' Assemblies are unpolished. Requires mamba, the external drive
#' (`/Volumes/lan_images4/assembly/`), the cache from `mito-flye-round1.sh`, and
#' `/Volumes/lan_images4/` (sup nanopore reads; fastp-trimmed MyKS-19-B
#' reads, streamed from the uncompressed tar; the 2026-10-19 output used
#' the raw Illumina reads). Takes ~2-3 h. Run from the repo root.
#' Reads default to the Dorado sup basecalls on the external drive; outputs
#' saved before 2026-10-23 used the Guppy 5.0.11 hac reads (deleted
#' locally), as do reads already in `_data/mito-flye-cache/reassemble/`, which are
#' reused if present (delete them to recruit sup reads).
#' Set ENV_DIR to reuse environments (`asm` and `smart` inside it).
#'
#' Output: `_data/mito-reassemble.txt`; reads and assemblies in
#' `_data/mito-flye-cache/reassemble/`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-reassemble.txt
CACHE=${ROOT}/_data/mito-flye-cache
RDIR=${CACHE}/reassemble
R1C=${CACHE}/with_opts_assembly.fasta
OLD=/Volumes/lan_images4/assembly
ONT_FASTQ=${ONT_FASTQ:-/Volumes/lan_images4/ont/basecalls_dorado-0.9.6_sup.fastq.gz}
ILL=/Volumes/lan_images4/ill/dna/trimmed/trimmed_MyKS-19-B_S18.tar
ILL_PAIRS=4000000
TARGET_BP=22500000  # ~500x of a 45-kb molecule
THREADS=$(sysctl -n hw.ncpu 2> /dev/null || nproc)
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/envs}
cd ${WORK}

[ -s ${R1C} ] || { echo "run _bash/tests/mito-flye-round1.sh first" 1>&2; exit 1; }
[ -d ${OLD} ] || { echo "external drive not mounted" 1>&2; exit 1; }
[ -s ${ILL} ] || { echo "missing ${ILL} (drive mounted?)" 1>&2; exit 1; }

# ---- environments
if [ ! -d ${ENV_DIR}/asm ]; then
    mamba create -y -q -p ${ENV_DIR}/asm -c conda-forge -c bioconda \
        flye=2.9.6 raven-assembler=1.8.3 nextdenovo=2.5.2 filtlong=0.3.1 \
        minimap2=2.31 samtools=1.24 seqtk=1.5 > /dev/null 2>&1 \
        && ${ENV_DIR}/asm/bin/python -m pip install -q paralleltask==0.1.1 \
        || { echo "could not install assemblers" 1>&2; exit 1; }
fi
if [ ! -d ${ENV_DIR}/smart ]; then
    # no osx-arm64 build
    CONDA_SUBDIR=osx-64 mamba create -y -q -p ${ENV_DIR}/smart \
        -c conda-forge -c bioconda smartdenovo=1.0.0 > /dev/null 2>&1 \
        || { echo "could not install smartdenovo" 1>&2; exit 1; }
fi
export PATH=${ENV_DIR}/asm/bin:${PATH}
mkdir -p ${RDIR}

# ---- reference pieces (from the round-1 contig)
. ${ROOT}/_bash/tests/mito-features.sh
locate_features ${R1C} || exit 1
samtools faidx r1c.fa $(region genes) | sed '1s/.*/>genes/' > genes.fa
samtools faidx r1c.fa $(region CR_left) | sed '1s/.*/>CR_left/' > pieces.fa
samtools faidx r1c.fa $(region CR_right) | sed '1s/.*/>CR_right/' >> pieces.fa
samtools faidx r1c.fa $(region unit) | sed '1s/.*/>unit/' > unit.fa
get_ctg () {  # get_ctg <fasta.gz> <contig> <new name>
    gunzip -c $1 | awk -v c=$2 '/^>/ {p = (substr($1, 2) == c)} p' \
        | sed "1s/.*/>$3/"
}
get_ctg ${OLD}/nextdenovo/contigs_next.fasta.gz ctg000860 next2019 > next2019.fa
get_ctg ${OLD}/smartdenovo/contigs_smart.fasta.gz utg319 smart2019 > smart2019.fa
cat r1c.fa next2019.fa smart2019.fa > bait.fa

# ---- 1. recruit (cached)
if [ ! -s ${RDIR}/mito_reads.fastq ]; then
    minimap2 -x map-ont --secondary=no -t ${THREADS} bait.fa ${ONT_FASTQ} \
        2> /dev/null | awk '$11 >= 1000 {print $1}' | sort -u > mito.names
    seqtk subseq ${ONT_FASTQ} mito.names > ${RDIR}/mito_reads.fastq
fi
MR=${RDIR}/mito_reads.fastq

# ---- 2. subsample
filtlong --target_bases ${TARGET_BP} --length_weight 10 ${MR} \
    2> /dev/null > sub.fq
seqtk seq -a sub.fq > sub.fa

# ---- 3. assemble
flye --nano-raw sub.fq --out-dir flye --threads ${THREADS} \
    --genome-size 45k --asm-coverage 200 > /dev/null 2>&1
echo "flye exit $?" > asm_status.txt
cp flye/assembly.fasta asm_flye.fa 2> /dev/null

raven -t ${THREADS} --disable-checkpoints sub.fq 2> /dev/null > asm_raven.fa
echo "raven exit $?" >> asm_status.txt

realpath sub.fq > input.fofn 2> /dev/null || echo "$(pwd)/sub.fq" > input.fofn
cat > nd.cfg << EOF
[General]
job_type = local
job_prefix = nd
task = all
rewrite = yes
deltmp = yes
parallel_jobs = 2
input_type = raw
read_type = ont
input_fofn = input.fofn
workdir = nd

[correct_option]
read_cutoff = 1k
genome_size = 45k
sort_options = -m 4g -t $(( THREADS / 2 ))
minimap2_options_raw = -t $(( THREADS / 2 ))
pa_correction = 2
correction_options = -p $(( THREADS / 2 ))

[assemble_option]
minimap2_options_cns = -t $(( THREADS / 2 ))
nextgraph_options = -a 1
EOF
nextDenovo nd.cfg > nd.log 2>&1
echo "nextdenovo exit $?" >> asm_status.txt
cp nd.log ${RDIR}/nextdenovo.log
cp nd/03.ctg_graph/nd.asm.fasta asm_nextdenovo.fa 2> /dev/null

mkdir smart && (cd smart \
    && PATH=${ENV_DIR}/smart/bin:${PATH} smartdenovo.pl -p sm -t ${THREADS} \
        -c 1 ../sub.fa > sm.mak 2> /dev/null \
    && PATH=${ENV_DIR}/smart/bin:${PATH} make -f sm.mak > ../smart.log 2>&1)
echo "smartdenovo exit $?" >> asm_status.txt
cp smart/sm.dmo.cns asm_smartdenovo.fa 2> /dev/null

for a in flye raven nextdenovo smartdenovo; do
    [ -s asm_${a}.fa ] && cp asm_${a}.fa ${RDIR}/
done

# ---- 4-5. per-contig layout and tests
# mito contigs: >= 5 kb of the gene region aligned
: > cands.fa
for a in flye raven nextdenovo smartdenovo; do
    [ -s asm_${a}.fa ] || continue
    minimap2 -c -x asm20 --secondary=yes -N 5 asm_${a}.fa genes.fa \
        2> /dev/null | awk '$4 - $3 >= 5000 {print $6}' | sort -u > ids.txt
    while read id; do
        awk -v c=${id} '/^>/ {p = (substr($1, 2) == c)} p' asm_${a}.fa \
            | sed "1s/.*/>${a}_${id}/" >> cands.fa
    done < ids.txt
done
cat next2019.fa >> cands.fa
samtools faidx cands.fa

# long-read and Illumina alignments to all candidates at once would make
# reads choose among near-identical candidates, so align per candidate
for k in 1 2; do
    tar -xOf ${ILL} trimmed_MyKS-19-B_S18_L002_R${k}_001.fastq.gz 2> /dev/null \
        | gunzip -c 2> /dev/null | head -n $(( ILL_PAIRS * 4 )) > ill${k}.fq
done

# span_counts <bam> <contig> <positions file>: reads spanning pos +/- 1 kb
span_counts () {
    samtools view -F 0x904 -q 20 $1 $2 \
        | awk 'NR == FNR {p[++n] = $1; next}
               {e = $4; c = $6
                while (match(c, /[0-9]+[MDN=X]/)) {
                    e += substr(c, RSTART, RLENGTH - 1); c = substr(c, RSTART + RLENGTH)
                    }
                for (i = 1; i <= n; i++) if ($4 <= p[i] - 1000 && e >= p[i] + 1000) s[i]++}
               END {for (i = 1; i <= n; i++) printf " %d:%d", p[i], s[i] + 0
                    printf "\n"}' $3 -
}

{
    echo "flye 2.9.6, raven 1.8.3, nextdenovo 2.5.2, smartdenovo 1.0.0," \
        "filtlong 0.3.1, minimap2 2.31, samtools 1.24, seqtk 1.5;" \
        "checked $(date +%Y-%m-%d)"
    echo "recruited reads: $(awk 'NR % 4 == 2 {n++; b += length($0)}
        END {printf "%d reads, %d bp", n, b}' ${MR}); subsample:" \
        "$(awk 'NR % 4 == 2 {n++; b += length($0); if (length($0) < m || !m)
        m = length($0)} END {printf "%d reads, %d bp, shortest %d", n, b, m}' sub.fq)"
    sed 's/^/assembler /' asm_status.txt
    for a in flye raven nextdenovo smartdenovo; do
        printf "%s contigs (name length):" ${a}
        awk '/^>/ {if (n) printf " %s %d", n, l; n = substr($1, 2); l = 0; next}
             {l += length($0)} END {if (n) printf " %s %d", n, l; printf "\n"}' \
            asm_${a}.fa 2> /dev/null
    done
    while read -u 3 cand len rest; do
        samtools faidx cands.fa ${cand} > cand.fa
        echo
        echo "== ${cand}: ${len} bp, A+T $(awk '!/^>/ {s = s $0} END {n = length(s)
            a = gsub(/[AaTt]/, "", s); printf "%.3f", a / n}' cand.fa)"
        echo "layout (piece start-end ~ contig start-end strand aligned-bp identity):"
        minimap2 -c -x asm20 --secondary=yes -N 20 -p 0.3 cand.fa genes.fa \
            pieces.fa 2> /dev/null \
            | awk '$11 >= 500 {printf "  %-8s %5d-%5d ~ %6d-%6d %s %5d %.3f\n",
                   $1, $3, $4, $8, $9, $5, $11, $10 / $11}' | sort -k3,3n
        minimap2 -c -x asm20 --secondary=yes -N 5 cand.fa genes.fa 2> /dev/null \
            | awk '$4 - $3 >= 15000 {print $8}' | sort -n \
            | awk 'NR > 1 {printf "gene-region starts %d apart\n", $1 - p} {p = $1}'
        minimap2 -c -x asm20 --secondary=yes -N 50 -p 0.3 cand.fa unit.fa \
            2> /dev/null | awk '$11 >= 1500 {print $8, $9}' | sort -n > units.txt
        echo "repeat-unit copies (>= 1.5 kb aligned): $(wc -l < units.txt | tr -d ' ')"
        # boundaries between neighbouring units, and between the unit array
        # and its flanks
        awk 'NR > 1 {print int((pe + $1) / 2)} {pe = $2}' units.txt > bounds.txt
        awk 'NR == 1 {print $1} END {print $2}' units.txt >> bounds.txt
        sort -n bounds.txt -o bounds.txt
        minimap2 -a -x map-ont --secondary=no -t ${THREADS} cand.fa ${MR} \
            2> /dev/null | samtools sort -o ont.bam - 2> /dev/null
        samtools index ont.bam
        if [ -s bounds.txt ]; then
            echo "reads spanning +/- 1 kb of each unit boundary (position:reads):"
            echo " $(span_counts ont.bam ${cand} bounds.txt)" | fold -s -w 76 \
                | sed 's/^/ /; s/[[:space:]]*$//'
        fi
        minimap2 -a -x sr -t ${THREADS} cand.fa ill1.fq ill2.fq 2> /dev/null \
            | samtools sort -o ill.bam - 2> /dev/null
        samtools index ill.bam
        samtools depth -a -Q 0 ill.bam > ill.depth
        GENE=$(minimap2 -c -x asm20 --secondary=no cand.fa genes.fa 2> /dev/null \
            | awk '{print $8, $9; exit}')
        if [ -n "${GENE}" ] && [ -s units.txt ]; then
            echo "Illumina depth (all MAPQ) per unit / gene-region depth; sum:"
            awk -v g="${GENE}" 'BEGIN {split(g, x, " ")}
                NR == FNR {s[NR] = $1; e[NR] = $2; n = NR; next}
                $2 > x[1] && $2 <= x[2] {gs += $3; gn++}
                {for (i = 1; i <= n; i++) if ($2 > s[i] && $2 <= e[i]) {us[i] += $3; un[i]++}}
                END {gd = gs / gn; printf "  gene region %.0fx;", gd
                     for (i = 1; i <= n; i++) {r = (us[i] / un[i]) / gd; t += r
                                               printf " %.2f", r}
                     printf "; sum %.1f\n", t}' units.txt ill.depth | fold -s -w 76 \
                | sed 's/[[:space:]]*$//'
        fi
    done 3< cands.fa.fai
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
