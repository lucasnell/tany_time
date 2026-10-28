#!/bin/bash

#' Build the reference for competitive mapping ("hologenome"; as in DEST,
#' Kapun et al. 2021): the T. gracilentus assembly plus assemblies of other
#' species whose DNA could be in the samples. Reads then map to whichever
#' genome they match best, and only reads assigned to T. gracilentus
#' contigs are used downstream (`03-mpileup`).
#'
#' T. gracilentus contig names are unchanged, except that the assembly's own
#' mitochondrial contigs (listed in `mt_contigs.txt` of the `01-mitogenome`
#' output: contigs carrying COI MZ659759.1; CONTIG34, a 49-kb contig with
#' the gene region twice) are dropped, so mitochondrial reads map only to
#' `mt_Tgraci__` and stay out of `target_contigs.bed` (owner's decision
#' 2026-10-28). Other sequences are renamed
#' so they're easy to filter out and count by group:
#' - `nt_<code>__<name>`: nuclear assemblies of other species
#' - `mt_<code>__<name>`: their mitogenomes, so their mitochondrial reads
#'   are counted apart from their nuclear reads (2026-10-28): the ones in
#'   the C. tentans and P. akamusi assemblies (OTHER_MT, relabeled), and
#'   NCBI mitogenomes for C. riparius and S. aterrima, whose assemblies
#'   have none (`other-mitogenomes.fasta` from `_R/get-other-mitogenomes.R`;
#'   `_bash/tests/other-mito-check.sh`)
#' - `mt_Tgraci__Tgraci_mitogenome`: T. gracilentus mitogenome
#'   (`01-mitogenome`), so mitochondrial reads don't map to nuclear copies
#'   of mitochondrial DNA
#' - `bc_<code>__<accession>`: COI barcodes of abundant trap taxa
#'   (`coi-barcodes.fasta` from `_R/get-mito-refs.R`). Reads on these vs.
#'   the T. gracilentus mitogenome give a per-pool contamination index.
#' Only `target_contigs.bed` (nuclear T. gracilentus) is used for mpileup.
#'
#' Output `hologenome.tar` (on staging) contains directory `hologenome/`:
#' - hologenome.fasta (+ bwa index, .fai, .dict)
#' - target_contigs.bed: T. gracilentus contigs (0-based, full length)
#'


. /app/.bashrc
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export GENOME_FULL_PATH="/staging/l/lnell/Tgraci_assembly.fasta.gz"
export OTHER_DIR="/staging/l/lnell/other-genomes"
export PARENT_DIR_OUT="/staging/l/lnell"
export MITO_TAR_FULL_PATH="/staging/l/lnell/Tgraci_mitogenome.tar"

#' Non-target assemblies, as `<code>:<file in OTHER_DIR>`, chosen as the
#' closest available genomes to abundant non-target taxa in Mývatn traps
#' (Einarsson et al. 2026, Aquatic Ecology 60:58, Table 1):
#' - Chironomini (Chironomus islandicus): Chironomus riparius
#'   (GCA_917627325) and C. tentans (GCA_963573255.1, 2023)
#' - Orthocladiinae (e.g., Cricotopus sylvestris, Eukiefferiella minor):
#'   Smittia aterrima (GCA_033063855.1) and Propsilocerus akamusi
#'   (GCA_018397935)
#' No assemblies exist for Micropsectra (Tanytarsini), Tanypodinae, or
#' Simulium (checked 2026-10-12).
export OTHERS=("Cripar:Cripar_assembly.fasta.gz"
               "Ctenta:Ctenta2023_assembly.fasta.gz"
               "Saterr:Saterr_assembly.fasta.gz"
               "Pakamu:Pakamu_assembly.fasta.gz")
#' Mitogenomes inside those assemblies, as `<code>:<sequence name>`
#' (`_bash/tests/other-mito-check.sh`):
export OTHER_MT=("Ctenta:OY754306.1" "Pakamu:JAFFZY010000103.1")


OTHER_FILES=()
for x in "${OTHERS[@]}"; do OTHER_FILES+=(${OTHER_DIR}/${x#*:}); done
require_files ${GENOME_FULL_PATH} ${OTHER_FILES[@]} ${MITO_TAR_FULL_PATH} \
    coi-barcodes.fasta other-mitogenomes.fasta


export OUT_DIR=hologenome
export HOLO=hologenome.fasta
mkdir ${OUT_DIR}
mv coi-barcodes.fasta other-mitogenomes.fasta ${OUT_DIR}/
cd ${OUT_DIR}


#' Target genome first, so its contigs come first in the sequence
#' dictionary. Headers are trimmed to their first word; mitochondrial
#' contigs are dropped.
tar -xf ${MITO_TAR_FULL_PATH} -O Tgraci_mitogenome/mt_contigs.txt > mt_contigs.txt
check_exit_status "read mt_contigs.txt" $?
echo "Dropping the assembly's mitochondrial contigs: $(paste -sd ' ' mt_contigs.txt)"
gunzip -c ${GENOME_FULL_PATH} \
    | awk 'NR == FNR {mt[$1]; next}
           /^>/ {sub(/^>/, ""); n = $1; keep = !(n in mt); if (keep) print ">" n; next}
           keep' mt_contigs.txt - \
    > ${HOLO}
check_exit_status "add target genome" $?
N_DROPPED=$(gunzip -c ${GENOME_FULL_PATH} | grep "^>" | sed 's/^>//; s/ .*//' \
    | grep -cxFf mt_contigs.txt)
echo "Dropped ${N_DROPPED} of $(wc -l < mt_contigs.txt) listed contigs"

samtools faidx ${HOLO}
check_exit_status "faidx (target)" $?
awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' ${HOLO}.fai > target_contigs.bed
rm ${HOLO}.fai

for x in "${OTHERS[@]}"; do
    code=${x%%:*}
    mt=$(for y in "${OTHER_MT[@]}"; do [ "${y%%:*}" != "${code}" ] || echo "${y#*:}"; done)
    gunzip -c ${OTHER_DIR}/${x#*:} \
        | awk -v code=${code} -v mt="${mt}" '
            BEGIN {n = split(mt, m, " "); for (i = 1; i <= n; i++) M[m[i]]}
            /^>/ {sub(/^>/, ""); print ">" (($1 in M) ? "mt_" : "nt_") code "__" $1; next}
            {print}' \
        >> ${HOLO}
    check_exit_status "add ${code}" $?
done

# NCBI mitogenomes for the other species without one in their assembly
# (already named `mt_<code>__<accession>`):
cat other-mitogenomes.fasta >> ${HOLO} && rm other-mitogenomes.fasta
check_exit_status "add other mitogenomes" $?

# T. gracilentus mitogenome:
tar -xf ${MITO_TAR_FULL_PATH} -O Tgraci_mitogenome/Tgraci_mitogenome.fasta \
    | awk '{if (/^>/) {print ">mt_Tgraci__Tgraci_mitogenome"} else {print}}' \
    >> ${HOLO}
check_exit_status "add mitogenome" $?

# COI barcodes (already named `bc_<code>__<accession>`):
BC_GROUPS=($(grep "^>" coi-barcodes.fasta | sed 's/^>//; s/__.*//' | sort -u))
cat coi-barcodes.fasta >> ${HOLO} && rm coi-barcodes.fasta
check_exit_status "add barcodes" $?

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

# Summary of sequences per group:
awk '{split($1, a, "__"); sp = (a[1] ~ /^(nt|mt|bc)_/) ? a[1] : "target";
     n[sp]++; len[sp] += $2}
     END {for (s in n) print s "\t" n[s] "\t" len[s]}' \
    ${HOLO}.fai | sort > group-summary.tsv
cat group-summary.tsv

# Check that every expected group is present (fail without output if not):
EXPECTED=("target" "mt_Tgraci" "${BC_GROUPS[@]}")
for x in "${OTHERS[@]}"; do EXPECTED+=("nt_${x%%:*}" "mt_${x%%:*}"); done
MISSING=0
for g in "${EXPECTED[@]}"; do
    if ! awk -v g="${g}" '$1 == g && $2 > 0 {f = 1} END {exit !f}' group-summary.tsv; then
        echo "MISSING from combined reference: ${g}" 1>&2
        MISSING=1
    fi
done
if (( MISSING != 0 )); then
    echo "HOLOGENOME CHECK FAILED. No output written." 1>&2
    safe_exit 1
fi
echo "HOLOGENOME CHECK PASSED: all ${#EXPECTED[@]} groups present."


cd ..
tar -cf ${OUT_DIR}.tar ${OUT_DIR}
check_exit_status "tar" $?
mv ${OUT_DIR}.tar ${PARENT_DIR_OUT}/
rm -r ${OUT_DIR}

safe_exit 0
