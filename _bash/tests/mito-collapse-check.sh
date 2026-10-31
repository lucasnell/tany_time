#!/bin/bash

#' Local check (not a cluster job): does `mito-collapse.sh`
#' (`bin/`, used by `00-references/01-mitogenome`) give the
#' expected layout on the mitochondrial contigs at hand?
#'
#' Contigs: the polished contig from the cluster run of `01-mitogenome`
#' (`Tgraci_mito_contig.fasta` in `_data/Tgraci_mitogenome.tar`, copied from
#' staging), the cached round-1 contigs (sup and Guppy), the Flye and
#' SMARTdenovo contigs from `mito-reassemble.sh`, and the 2019 NextDenovo
#' and SMARTdenovo contigs. For each, reports the layout, the output length
#' against the sum of the pieces, how much of the `genes` feature the output
#' carries (minimap2 -c -x asm20), and the COI barcodes' identity. Contigs
#' without a full gene region or a full unit should fail.
#' Needs minimap2 and samtools: set ENV_DIR to an environment with both.
#' Requires `/Volumes/lan_images4/assembly/`. Run from the repo root.
#'
#' Output: `_data/mito-collapse-check.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/mito-collapse-check.txt
CACHE=${ROOT}/_data/mito-flye-cache
FEATS=${ROOT}/_bash/references/tgraci-mito-features.fasta
COI=${ROOT}/_bash/references/tgraci-coi.fasta
COLLAPSE=${ROOT}/_bash/bin/mito-collapse.sh
OLD=/Volumes/lan_images4/assembly
MT_TAR=${ROOT}/_data/Tgraci_mitogenome.tar
WORK=$(mktemp -d)
[ -n "${ENV_DIR}" ] && export PATH=${ENV_DIR}/bin:${PATH}
for t in minimap2 samtools; do
    command -v ${t} > /dev/null || { echo "${t} not found (set ENV_DIR)" 1>&2; exit 1; }
done
[ -d ${OLD} ] || { echo "missing ${OLD} (drive mounted?)" 1>&2; exit 1; }
cd ${WORK}
cp ${FEATS} feats.fa

# one sequence per file: label, source
get () {  # get <fasta(.gz)> <contig> <out>
    gunzip -cf $1 | awk -v c=$2 '/^>/ {p = (substr($1, 2) == c)} p' > $3
}
tar -xOf ${MT_TAR} Tgraci_mitogenome/Tgraci_mito_contig.fasta > cluster.fa
get ${CACHE}/with_opts_assembly.fasta contig_1 r1_sup.fa
get ${CACHE}-guppy/with_opts_assembly.fasta contig_10 r1_guppy.fa
get ${CACHE}/reassemble/asm_flye.fa contig_2 flye2.fa
get ${CACHE}/reassemble/asm_flye.fa contig_3 flye3.fa
get ${CACHE}/reassemble/asm_smartdenovo.fa utg363 smart_sup.fa
get ${OLD}/nextdenovo/contigs_next.fasta.gz ctg000860 next2019.fa
get ${OLD}/smartdenovo/contigs_smart.fasta.gz utg319 smart2019.fa

{
    echo "minimap2 $(minimap2 --version); checked $(date +%Y-%m-%d)"
    for x in cluster r1_sup r1_guppy flye2 flye3 smart_sup next2019 smart2019; do
        echo
        echo "== ${x}: $(head -1 ${x}.fa | sed 's/^>//; s/ .*//')," \
            "$(awk '!/^>/ {n += length($0)} END {print n}' ${x}.fa) bp"
        if ! bash ${COLLAPSE} ${x}.fa feats.fa ${x}_mt 2>&1; then
            echo "  collapse failed"; continue
        fi
        sed 's/^/  /' ${x}_mt_layout.tsv | column -t | sed 's/^/  /'
        n=$(awk '!/^>/ {n += length($0)} END {print n}' ${x}_mt.fasta)
        p=$(awk 'NR > 1 {s += $6} END {print s}' ${x}_mt_layout.tsv)
        echo "  output ${n} bp; pieces sum ${p} bp: $([ ${n} -eq ${p} ] && echo OK || echo MISMATCH)"
        samtools faidx feats.fa genes \
            | minimap2 -c -x asm20 ${x}_mt.fasta - 2> /dev/null \
            | sort -k10,10nr | awk 'NR == 1 {printf "  genes feature: %d of %d bp aligned, identity %.4f\n",
                $4 - $3, $2, $10 / $11}'
        minimap2 -c -x asm20 --secondary=no ${x}_mt.fasta ${COI} 2> /dev/null \
            | awk '{printf "  COI %s: %d of %d bp, identity %.4f\n", $1, $4 - $3, $2, $10 / $11}'
    done
} > ${OUT}
sed -i '' 's/[[:space:]]*$//' ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
