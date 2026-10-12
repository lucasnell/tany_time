#!/bin/bash

#' Build the reference for competitive mapping ("hologenome"; as in DEST,
#' Kapun et al. 2021): the T. gracilentus assembly plus assemblies of other
#' species whose DNA could be in the samples. Reads then map to whichever
#' genome they match best, and only reads assigned to T. gracilentus
#' contigs are used downstream (`03-mpileup`).
#'
#' T. gracilentus contig names are unchanged. Other sequences are renamed
#' `nt_<code>__<first word of header>` so they're easy to filter out.
#'
#' Output `hologenome.tar` (on staging) contains directory `hologenome/`:
#' - hologenome.fasta (+ bwa index, .fai, .dict)
#' - target_contigs.bed: T. gracilentus contigs (0-based, full length)
#'


. /app/.bashrc
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export GENOME_FULL_PATH="/staging/lnell/Tgraci_assembly.fasta.gz"
export OTHER_DIR="/staging/lnell/other-genomes"
export PARENT_DIR_OUT="/staging/lnell"

#' Non-target assemblies, as `<code>:<file in OTHER_DIR>`. Chironomus spp.
#' are proxies for the Chironomus species at Mývatn (no assembly available).
export OTHERS=("Cripar:Cripar_assembly.fasta.gz"
               "Ctenta:Ctenta_assembly.fasta.gz")


OTHER_FILES=()
for x in "${OTHERS[@]}"; do OTHER_FILES+=(${OTHER_DIR}/${x#*:}); done
require_files ${GENOME_FULL_PATH} ${OTHER_FILES[@]}


export OUT_DIR=hologenome
export HOLO=hologenome.fasta
mkdir ${OUT_DIR}
cd ${OUT_DIR}


#' Target genome first, so its contigs come first in the sequence
#' dictionary. Headers are trimmed to their first word.
gunzip -c ${GENOME_FULL_PATH} \
    | awk '{if (/^>/) {print $1} else {print}}' \
    > ${HOLO}
check_exit_status "add target genome" $?

samtools faidx ${HOLO}
check_exit_status "faidx (target)" $?
awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' ${HOLO}.fai > target_contigs.bed
rm ${HOLO}.fai

for x in "${OTHERS[@]}"; do
    code=${x%%:*}
    gunzip -c ${OTHER_DIR}/${x#*:} \
        | awk -v pre="nt_${code}__" \
            '{if (/^>/) {sub(/^>/, ""); print ">" pre $1} else {print}}' \
        >> ${HOLO}
    check_exit_status "add ${code}" $?
done

# Check that names are unique:
N_DUPS=$(grep "^>" ${HOLO} | sort | uniq -d | wc -l)
if (( N_DUPS > 0 )); then
    echo "Duplicate sequence names in ${HOLO}" 1>&2
    safe_exit 1
fi

samtools faidx ${HOLO}
check_exit_status "faidx" $?

bwa index ${HOLO}
check_exit_status "bwa index" $?

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
conda activate gatk-env
gatk CreateSequenceDictionary R=${HOLO}
check_exit_status "CreateSequenceDictionary" $?
conda deactivate
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Summary of sequences per species:
awk '{split($1, a, "__"); sp = (a[1] ~ /^nt_/) ? substr(a[1], 4) : "target";
     n[sp]++; len[sp] += $2}
     END {for (s in n) print s "\t" n[s] " sequences\t" len[s] " bp"}' \
    ${HOLO}.fai


cd ..
tar -cf ${OUT_DIR}.tar ${OUT_DIR}
check_exit_status "tar" $?
mv ${OUT_DIR}.tar ${PARENT_DIR_OUT}/
rm -r ${OUT_DIR}

safe_exit 0
