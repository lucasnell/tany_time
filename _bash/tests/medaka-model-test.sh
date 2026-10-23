#!/bin/bash

#' Local check (not a cluster job): does medaka `r941_min_sup_g507` help or
#' hurt when polishing with the Dorado 0.9.6 sup v3.6 reads? Medaka has no
#' model for `dna_r9.4.1_e8_sup@v3.6`; its table pairs `r941_min_sup_g507`
#' with `dna_r9.4.1_e8_sup@v3.3` (`medaka/options.py`).
#'
#' Polishes the cached round-1 contig (sup reads; `mito-flye-round1.sh`)
#' with the same recruited reads and each model in MODELS, and compares
#' those and the unpolished draft against Illumina DNA reads (first 4M
#' trimmed pairs of pool MyKS-19-B_S18, cut from
#' `/Volumes/lan_images4/ill/dna/trimmed/trimmed_MyKS-19-B_S18.tar` into
#' `_data/mito-flye-cache/ill/`). For each version, Illumina reads are
#' mapped to it (minimap2 -x sr) and called against it (bcftools, haploid;
#' MAPQ >= 20, base quality >= 20; `-L` raised so indels are called at
#' mitochondrial depth, above bcftools' default 250). At sites with depth >= MIN_DP, a call
#' whose alternative allele has >= MIN_AF of the reads is an error in that
#' version; MIX_AF to MIN_AF is counted as mixed (pool polymorphism or
#' ambiguous). Indels are homopolymer indels if the inserted or deleted
#' bases are one base repeated and extend a run of >= 3 of it. Rates are
#' per 10 kb of sites with depth >= MIN_DP, in the gene region (located by
#' `mito-features.sh`) and outside it (control-region side; the A+T-rich
#' repeat copies are mostly below MIN_DP).
#' The Illumina pool is not the nanopore individual, so true differences
#' between that individual's haplotype and the pool's majority count as
#' errors in every version equally; compare versions, not absolute rates.
#' Versions are read from `_apptainer/mito-env.yml` (bcftools from
#' `main-env.yml`). Requires mamba and internet (medaka models). Run from
#' the repo root. Set ENV_DIR to reuse an environment.
#'
#' Output: `_data/medaka-model-test.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/medaka-model-test.txt
YML=${ROOT}/_bash/_apptainer/mito-env.yml
YML_MAIN=${ROOT}/_bash/_apptainer/main-env.yml
CACHE=${ROOT}/_data/mito-flye-cache
R1C=${CACHE}/with_opts_assembly.fasta
READS=${CACHE}/reads1.fastq
ILL1=${CACHE}/ill/MyKS-19-B_S18_R1.fq
ILL2=${CACHE}/ill/MyKS-19-B_S18_R2.fq
MODELS=(r941_min_hac_g507 r941_min_sup_g507)
MIN_DP=50
MIN_AF=0.9
MIX_AF=0.1
THREADS=$(sysctl -n hw.ncpu 2> /dev/null || nproc)
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/env}
cd ${WORK}

for f in ${R1C} ${READS} ${ILL1} ${ILL2}; do
    [ -s ${f} ] || { echo "missing ${f}" 1>&2; exit 1; }
done

pin () { grep -E "^\s+- $1=" $2 | sed 's/.*=//'; }
PKGS="medaka=$(pin medaka ${YML}) minimap2=$(pin minimap2 ${YML})
      samtools=$(pin samtools ${YML}) bcftools=$(pin bcftools ${YML_MAIN})"

if [ ! -x ${ENV_DIR}/bin/medaka ]; then
    mamba create -y -q -p ${ENV_DIR} -c conda-forge -c bioconda ${PKGS} \
        > /dev/null 2>&1 \
        || CONDA_SUBDIR=osx-64 mamba create -y -q -p ${ENV_DIR} \
            -c conda-forge -c bioconda ${PKGS} > /dev/null 2>&1 \
        || { echo "could not install ${PKGS}" 1>&2; exit 1; }
fi
export PATH=${ENV_DIR}/bin:${PATH}
medaka tools download_models --models ${MODELS[@]} > /dev/null 2>&1 \
    || { echo "could not download medaka models" 1>&2; exit 1; }

. ${ROOT}/_bash/tests/mito-features.sh
locate_features ${R1C} || exit 1
mv r1c.fa draft.fa

VERSIONS=(draft)
for m in ${MODELS[@]}; do
    medaka_consensus -i ${READS} -d draft.fa -o med_${m} -m ${m} \
        -t ${THREADS} > med_${m}.log 2>&1 \
        || { echo "medaka ${m} failed; see ${WORK}/med_${m}.log" 1>&2; exit 1; }
    awk '/^>/ && n++ {exit} {print}' med_${m}/consensus.fasta > ${m}.fa
    VERSIONS+=(${m})
done

# score <version>: one line of the results table
score () {
    local v=$1
    mkdir -p s_${v} && cd s_${v}
    locate_features ../${v}.fa || return 1
    samtools faidx r1c.fa
    minimap2 -a -x sr -t ${THREADS} r1c.fa ${ILL1} ${ILL2} 2> /dev/null \
        | samtools sort -@ 2 -o ill.bam - 2> /dev/null
    samtools index ill.bam
    samtools depth -a -Q 20 -q 20 ill.bam \
        | awk -v m=${MIN_DP} -v gs=${GS} -v ge=${GE} '$3 >= m {
              if ($2 >= gs && $2 <= ge) g++; else o++}
              END {print g + 0, o + 0}' > callable.txt
    bcftools mpileup -Ou -f r1c.fa -d 100000 -L 100000 -q 20 -Q 20 -a AD,DP ill.bam \
        2> /dev/null | bcftools call -m -v --ploidy 1 -Ov 2> /dev/null \
        | bcftools query -f '%POS\t%REF\t%ALT\t%INFO/DP\t[%AD]\n' > calls.txt
    # sequence on one line for homopolymer checks
    awk '!/^>/ {printf "%s", $0} END {print ""}' r1c.fa > seq.txt
    # repeat-unit copies (hits of the unit covering >= 1.5 kb)
    samtools faidx ${ROOT}/_data/mito-features.fasta unit \
        | minimap2 -c -x asm20 --secondary=yes -N 50 -p 0.3 r1c.fa - 2> /dev/null \
        | awk '$4 - $3 >= 1500 {printf "%s%d-%d", (n++ ? "," : ""), $8 + 1, $9}' > units.txt
    awk -v m=${MIN_DP} -v u=$(cat units.txt) -v v=${v} -v hi=${MIN_AF} -v lo=${MIX_AF} -v gs=${GS} -v ge=${GE} \
        -v c=$(cat callable.txt | tr ' ' ,) '
        BEGIN {nu = split(u, U, ",")}
        NR == FNR {s = toupper($0); next}
        {split($5, ad, ","); t = 0; for (i in ad) t += ad[i]
         if (t < m || $4 < m) next
         alt = ad[2] / t; r = ($1 >= gs && $1 <= ge) ? "g" : "o"
         if (alt < lo) next
         if (alt < hi) {mix[r]++; next}
         split($3, A, ","); a = A[1]
         inu = "no"; for (i = 1; i <= nu; i++) {split(U[i], ue, "-")
             if ($1 >= ue[1] && $1 <= ue[2]) inu = "yes"}
         site = sprintf("  %-18s %-7s %6d %s>%s alt %.2f depth %d in-unit %s",
             v, (r == "g" ? "genes" : "outside"), $1, $2, a, alt, t, inu)
         if (length($2) == length(a)) {snv[r]++; print site " SNV" >> "../sites.txt"; next}
         # inserted or deleted bases, after the shared first base
         d = (length($2) > length(a)) ? substr($2, 2) : substr(a, 2)
         b = substr(d, 1, 1); hp = (d ~ ("^" b "+$"))
         # run of b in the reference around POS (deleted bases included;
         # plus the inserted bases for insertions)
         k = 0; while (substr(s, $1 + 1 + k, 1) == b) k++
         j = 0; while (substr(s, $1 - j, 1) == b) j++
         if (hp && (k + j + (length(a) > length($2) ? length(d) : 0)) >= 3) {
             hpi[r]++; print site " homopolymer-indel" >> "../sites.txt"}
         else {oi[r]++; print site " other-indel" >> "../sites.txt"}}
        END {split(c, cc, ",")
             for (j = 1; j <= 2; j++) {r = (j == 1) ? "g" : "o"
                 e = snv[r] + hpi[r] + oi[r]
                 printf " | %6d %4d %4d %4d %4d %6.2f", cc[j], snv[r], hpi[r],
                     oi[r], mix[r], (cc[j] ? 1e4 * e / cc[j] : 0)}}' \
        seq.txt calls.txt > row.txt
    printf "%-18s %6d%s\n" ${v} ${L} "$(cat row.txt)"
    cd ..
}

{
    echo "medaka $(medaka --version 2>&1 | awk '{print $NF}'), minimap2" \
        "$(minimap2 --version), samtools $(pin samtools ${YML})," \
        "bcftools $(pin bcftools ${YML_MAIN}); checked $(date +%Y-%m-%d)"
    echo "draft: ${R1C} (${L} bp; gene region $(region genes));" \
        "reads: $(awk 'END {print NR / 4}' ${READS}) recruited sup reads"
    echo "Illumina: $(awk 'END {print NR / 4}' ${ILL1}) pairs; sites with" \
        "depth >= ${MIN_DP}; error = alt >= ${MIN_AF} of reads; mixed =" \
        "${MIX_AF}-${MIN_AF}"
    echo
    echo "columns per region: callable bp, SNV, homopolymer indel, other" \
        "indel, mixed, errors per 10 kb"
    printf "%-18s %6s | %-39s | %s\n" version length \
        "gene region" "outside gene region"
    : > sites.txt
    for v in ${VERSIONS[@]}; do
        score ${v} || { echo "scoring ${v} failed" 1>&2; exit 1; }
    done
    echo
    echo "error sites (version, region, position on that version, ref>alt," \
        "alt fraction, depth, inside a repeat-unit copy, type):"
    cat sites.txt
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
