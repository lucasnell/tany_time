#!/bin/bash

#' Check a pipeline step after (or while) it runs, on the submit server.
#'
#' Usage (from anywhere): `_bash/check-step.sh <folder>/<step>`, e.g.,
#'   _bash/check-step.sh 01-variant-calling/02-bwa
#'
#' Reports:
#' 1. Jobs: finished with exit code 0, failed, still running or held, and
#'    expected jobs with no HTCondor log (not submitted).
#' 2. Run time: longest job (or longest so far), with a warning above
#'    WARN_H hours (CHTC's default limit is 72 h).
#' 3. Expected outputs missing from staging.
#' 4. A step-specific summary (e.g., COI check, reads per species, SNPs
#'    kept by BayPass input filters, Ω consistency).
#'
#' Ends with `RESULT: PASS` (exit 0) or `RESULT: <n> problem(s)` (exit 1).
#' Warnings don't count as problems. Uses GNU `date` (or python3) to read
#' HTCondor log timestamps.
#'

STEP_PATH=${1%/}
if [ -z "${STEP_PATH}" ]; then
    echo "Usage: $0 <folder>/<step>, e.g., $0 01-variant-calling/02-bwa" 1>&2
    exit 2
fi

# Staging folder (override with the STAGING environment variable):
S=${STAGING:-/staging/l/lnell}
WARN_H=54

BASH_DIR=$(cd "$(dirname "$0")" && pwd)
FOLDER=${STEP_PATH%%/*}
STEP=${STEP_PATH##*/}
LOG_DIR=${BASH_DIR}/${FOLDER}/logs/${STEP}
SAMPLES=${BASH_DIR}/samples
ST=${BASH_DIR}/03-select-time

PROBLEMS=0
problem () { echo "PROBLEM: $*"; PROBLEMS=$((PROBLEMS + 1)); }
warn () { echo "WARN: $*"; }

if [ ! -d "${LOG_DIR}" ]; then
    echo "No log folder ${LOG_DIR} (submit from inside ${FOLDER}/)" 1>&2
    exit 2
fi

pools () { cat ${SAMPLES}/all-pools.txt; }

# Seconds since epoch for an HTCondor log event line
# ("NNN (cluster.proc.sub) YYYY-MM-DD HH:MM:SS ...").
event_secs () {
    local ts=$(echo "$1" | awk '{print $3, $4}')
    date -d "${ts}" +%s 2> /dev/null \
        || python3 -c "import sys, datetime; print(int(datetime.datetime.strptime(sys.argv[1], '%Y-%m-%d %H:%M:%S').timestamp()))" "${ts}" 2> /dev/null
}



# =============================================================================
# Step definitions: job IDs (HTCondor log names), expected outputs, extras
# =============================================================================

JOB_IDS=()
FILES=()
extra () { :; }

case "${STEP_PATH}" in

    00-basecall/01-split-fast5)
        JOB_IDS=(split-fast5)
        FILES=(${S}/dna/ont/fast5_chunks/manifest.tsv)
        extra () {
            grep -h "fast5 files:\|archive MD5:" ${LOG_DIR}/split-fast5.out 2> /dev/null
            local m=$(grep -h "archive MD5:" ${LOG_DIR}/split-fast5.out 2> /dev/null)
            [ -n "${m}" ] && echo "${m}" | awk '{exit !($3 == $5)}' \
                || problem "archive MD5 missing or not as expected"
            local mf=${S}/dna/ont/fast5_chunks/manifest.tsv
            if [ -f ${mf} ]; then
                local n_chunk=$(cut -f 1 ${mf} | sort -u | wc -l | tr -d " ")
                local n_tar=$(ls ${S}/dna/ont/fast5_chunks/chunk_*.tar 2> /dev/null | wc -l | tr -d " ")
                echo "chunks in manifest: ${n_chunk}; chunk tars on staging: ${n_tar}"
                (( n_chunk == n_tar )) || problem "manifest and chunk tars disagree"
            fi
        } ;;

    00-basecall/02-dorado)
        for f in ${S}/dna/ont/fast5_chunks/chunk_*.tar; do
            [ -f "${f}" ] || continue
            c=$(basename ${f} .tar)
            JOB_IDS+=(${c})
            FILES+=(${S}/dna/ont/dorado/${c}.bam ${S}/dna/ont/dorado/${c}_summary.tsv)
        done
        extra () {
            local n=0 r
            for id in "${JOB_IDS[@]}"; do
                r=$(grep -h "Simplex reads basecalled" ${LOG_DIR}/${id}.out 2> /dev/null \
                    | grep -oE "[0-9]+$" | tail -n 1)
                n=$((n + ${r:-0}))
            done
            echo "reads basecalled (all chunks with logs): ${n}"
            grep -h "Samples/s" ${LOG_DIR}/*.out 2> /dev/null \
                | grep -oE "[0-9.]+e\+[0-9]+" \
                | awk '{s += $1; k++} END {if (k) printf "mean speed: %.3g samples/s over %d jobs\n", s / k, k}'
        } ;;

    00-basecall/03-combine)
        JOB_IDS=(combine)
        FILES=(${S}/dna/ont/basecalls_dorado-0.9.6_sup.fastq.gz
               ${S}/dna/ont/basecalls_dorado-0.9.6_sup--summary.tsv.gz)
        extra () {
            grep -h "chunks:\|all reads:\|FASTQ reads:" ${LOG_DIR}/combine.out 2> /dev/null \
                || problem "no totals in combine.out"
        } ;;

    00-references/01-mitogenome)
        JOB_IDS=(mitogenome)
        FILES=(${S}/Tgraci_mitogenome.tar)
        extra () {
            echo "collapsed layout (piece, contig start/end, output start/end, length, truncated):"
            awk '/^Collapsed layout/ {p = 1; next} p && /^(piece|genes|unique|unit)/ {print} p && !/^(piece|genes|unique|unit)/ {p = 0}' \
                ${LOG_DIR}/mitogenome.out 2> /dev/null
            grep -q "^unit" ${LOG_DIR}/mitogenome.out 2> /dev/null \
                && grep -qE "^unique_(left|right).*yes$" ${LOG_DIR}/mitogenome.out 2> /dev/null \
                && warn "a unique control-region side reached the contig's end (truncated)"
            echo "COI check (barcode, length, aligned bp, identity, mitogenome length):"
            grep -hE "^(barcode|M[NZ][0-9]+)" ${LOG_DIR}/mitogenome.out 2> /dev/null
            grep -h "COI CHECK" ${LOG_DIR}/mitogenome.{out,err} 2> /dev/null
            grep -q "COI CHECK PASSED" ${LOG_DIR}/mitogenome.out 2> /dev/null \
                || problem "COI check did not pass"
        } ;;

    00-references/02-hologenome)
        JOB_IDS=(hologenome)
        FILES=(${S}/hologenome.tar)
        extra () {
            grep -h "^Dropp" ${LOG_DIR}/hologenome.out 2> /dev/null
            echo "Groups (group, sequences, bp):"
            awk -F'\t' 'NF == 3 && $1 ~ /^(target|nt_|mt_|bc_)/' \
                ${LOG_DIR}/hologenome.out 2> /dev/null
            grep -h "HOLOGENOME CHECK\|MISSING" ${LOG_DIR}/hologenome.{out,err} 2> /dev/null
            grep -q "HOLOGENOME CHECK PASSED" ${LOG_DIR}/hologenome.out 2> /dev/null \
                || problem "combined-reference check did not pass"
        } ;;

    00-references/03-silva-db)
        JOB_IDS=(silva-db)
        FILES=(${S}/kraken2-db/silva-138.2.tar.gz) ;;

    00-references/04-ont-pi)
        JOB_IDS=(ont-pi)
        BC=$(grep "^export BASECALLS=" ${BASH_DIR}/00-references/04-ont-pi.sh | cut -d= -f2)
        FILES=(${S}/dna/ont/ont_variants_${BC}.windowed.pi.gz ${S}/dna/ont/ont_variants_${BC}.vcf.gz)
        extra () {
            local theta=$(grep -E "^[0-9.]+(e-?[0-9]+)?$" ${LOG_DIR}/ont-pi.out 2> /dev/null | tail -n 1)
            local prior=$(grep "^export THETA=" ${BASH_DIR}/02-allele-freqs/01-snape.sh | cut -d= -f2)
            if [ -z "${theta}" ]; then
                problem "no θ estimate found in ont-pi.out"
            else
                local rel=$(awk -v t=${theta} -v p=${prior} 'BEGIN {printf "%.1f", 100 * (t - p) / p}')
                echo "θ from nanopore reads: ${theta}; SNAPE prior in 01-snape.sh: ${prior} (${rel}% difference)"
                awk -v r=${rel} 'BEGIN {exit !(r > 10 || r < -10)}' \
                    && warn "θ differs from the SNAPE prior by more than 10%; consider updating THETA in 02-allele-freqs/01-snape.sh"
            fi
        } ;;

    01-variant-calling/01-fastp)
        JOB_IDS=($(pools))
        for p in ${JOB_IDS[@]}; do
            FILES+=(${S}/dna/trimmed/trimmed_${p}.tar ${S}/dna/trimmed/trimmed_${p}.tar.gz)
        done ;;

    01-variant-calling/02-bwa)
        JOB_IDS=($(pools))
        for p in ${JOB_IDS[@]}; do
            FILES+=(${S}/dna/bwa/${p}_bwa.bam ${S}/dna/bwa/${p}_bwa.tar.gz)
        done
        extra () {
            echo "Reads per group (MAPQ >= 20, primary), as shares of mapped reads:"
            printf "%-16s %12s %8s %8s %8s %8s %8s\n" pool mapped target nt_other \
                mt_Tgraci mt_other barcode
            for p in ${JOB_IDS[@]}; do
                [ -f ${LOG_DIR}/${p}.out ] || continue
                awk -v p=${p} '
                    /^group\tlength\tn_reads/ {on = 1; next}
                    on && /^(target|nt_|mt_|bc_)[^\t]*\t[0-9.e+]+\t[0-9]+$/ {
                        g = ($1 ~ /^nt_/) ? "nt" : ($1 == "mt_Tgraci") ? "mt" : ($1 ~ /^mt_/) ? "mto" : ($1 ~ /^bc_/) ? "bc" : "target"
                        n[g] += $3
                        tot += $3; next }
                    on {on = 0}
                    END { if (tot > 0) printf "%-16s %12d %8.4f %8.4f %8.4f %8.5f %8.5f\n", p, tot,
                          n["target"] / tot, n["nt"] / tot, n["mt"] / tot, n["mto"] / tot, n["bc"] / tot }
                ' ${LOG_DIR}/${p}.out
            done | sort -k3,3n
            echo "(target includes mitochondrial reads on CONTIG37 in" \
                "references/mito-like-nuclear.bed; 03-pool-qc removes them)"
        } ;;

    01-variant-calling/03-mpileup)
        JOB_IDS=($(pools))
        for p in ${JOB_IDS[@]}; do
            FILES+=(${S}/dna/mpileup/${p}_mpileup.txt.gz ${S}/dna/mpileup/${p}_mpileup.tar.gz)
        done
        extra () {
            # bamtools stats prints "Total reads:" after duplicate removal,
            # then after indel realignment, which shouldn't change it.
            echo "Reads after duplicate removal and after realignment; pileup size:"
            printf "%-16s %12s %12s %10s\n" pool nodups realigned pileup_MB
            local p n1 n2 mb rows="" small="" msgs=()
            for p in ${JOB_IDS[@]}; do
                [ -f ${LOG_DIR}/${p}.out ] || continue
                read n1 n2 <<< $(awk '/^Total reads:/ {printf "%s ", $3}' ${LOG_DIR}/${p}.out)
                mb=NA
                [ -f ${S}/dna/mpileup/${p}_mpileup.txt.gz ] \
                    && mb=$(( $(stat -c %s ${S}/dna/mpileup/${p}_mpileup.txt.gz 2> /dev/null || stat -f %z ${S}/dna/mpileup/${p}_mpileup.txt.gz) / 1000000 ))
                rows+=$(printf "%-16s %12s %12s %10s" ${p} ${n1:-NA} ${n2:-NA} ${mb})$'\n'
                if [ -z "${n2}" ]; then
                    msgs+=("${p}: read counts missing from ${p}.out")
                elif [ "${n1}" != "${n2}" ] || [ "${n2}" = 0 ]; then
                    msgs+=("${p}: reads changed in realignment or are zero")
                fi
                [ "${mb}" = 0 ] && small+=" ${p}"
            done
            printf "%s" "${rows}" | sort -k3,3n
            for m in "${msgs[@]}"; do problem "${m}"; done
            # Near-empty pools (e.g., KS-2-83) are expected; pool QC drops them.
            [ -n "${small}" ] && warn "pileups under 1 MB:${small}"
        } ;;

    01-variant-calling/04-bam-cov)
        JOB_IDS=($(pools))
        for p in ${JOB_IDS[@]}; do FILES+=(${S}/dna/bwa/${p}_bwa_cov.bed.gz); done ;;

    01-variant-calling/05-kraken)
        JOB_IDS=(kraken)
        FILES=(${S}/dna/kraken/kraken.tar)
        extra () {
            local n_pools=$(pools | wc -l)
            local n_silva=$(grep -c "Checked step kraken2 SILVA" ${LOG_DIR}/kraken.out 2> /dev/null)
            local n_plus=$(grep -c "Checked step kraken2 (" ${LOG_DIR}/kraken.out 2> /dev/null)
            echo "Pools classified: SILVA ${n_silva:-0}, PlusPF ${n_plus:-0} (of ${n_pools})"
            grep -h "WARNING" ${LOG_DIR}/kraken.err 2> /dev/null
        } ;;

    02-allele-freqs/01-snape)
        JOB_IDS=($(pools))
        for p in ${JOB_IDS[@]}; do FILES+=(${S}/dna/snape/${p}_snape.tar); done
        extra () {
            local f=$(grep -l "Command failed for" ${LOG_DIR}/*.err 2> /dev/null)
            [ -n "${f}" ] && warn "contigs where SNAPE reported errors (see each pool's <contig>-*.err in its output tar):" \
                && grep -h "Command failed for" ${LOG_DIR}/*.err | sort | uniq -c
        } ;;

    02-allele-freqs/02-npstat)
        LOG_SUFFIX="_10kb"
        JOB_IDS=($(pools))
        for p in ${JOB_IDS[@]}; do FILES+=(${S}/dna/npstat/${p}_npstat_10kb.stat.gz); done
        extra () {
            local f=$(grep -l "Command failed for" ${LOG_DIR}/*.err 2> /dev/null)
            [ -n "${f}" ] && warn "contigs where npstat failed (1 = contig missing, 2 = npstat error):" \
                && grep -h "Command failed for" ${LOG_DIR}/*.err | sort | uniq -c
        } ;;

    02-allele-freqs/03-pool-qc)
        JOB_IDS=(pool-qc)
        FILES=(${S}/dna/pool-qc/pool-qc.tar)
        extra () {
            grep -h "WARNING" ${LOG_DIR}/pool-qc.err 2> /dev/null
            echo "Read shares (target_rate excludes mitolike_rate), lowest target first:"
            tar -xOf ${S}/dna/pool-qc/pool-qc.tar pool-qc/pool-qc-files.tsv 2> /dev/null \
                | awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i
                       printf "%-16s %9s %9s %9s\n", "pool", "map_rate", "target", "mitolike"; next}
                     {printf "%-16s %9s %9s %9s\n", $c["pool"], $c["map_rate"],
                         $c["target_rate"], $c["mitolike_rate"]}' \
                | { read -r h; echo "$h"; sort -k3,3g; }
        } ;;

    02-allele-freqs/04-combine-syncs)
        JOB_IDS=(combine-syncs)
        for k in snape count; do
            for t in time space all; do FILES+=(${S}/dna/snape/${k}_masked_${t}.tar); done
        done
        extra () {
            echo "Sites written (in order: snape time, space, all; count time, space, all):"
            grep -h "^wrote" ${LOG_DIR}/combine-syncs.out 2> /dev/null
        } ;;

    03-select-time/00-baypass-inputs)
        TAGS=($(cut -d, -f1 ${ST}/baypass_inputs_info.txt))
        for t in ${TAGS[@]}; do
            JOB_IDS+=(baypass-inputs_${t})
            FILES+=(${S}/dna/baypass/baypass-inputs_${t}.tar.gz ${S}/dna/baypass/baypass-inputs_${t}_3.tar.gz)
        done
        extra () {
            local prev=""
            for t in ${TAGS[@]}; do
                local f=${LOG_DIR}/baypass-inputs_${t}.out
                local kept=$(grep -h "^SNPs kept:" ${f} 2> /dev/null | head -n 1 | awk '{print $3}')
                echo "${t}: $(grep -h "^SNPs read:" ${f} 2> /dev/null | head -n 1); SNPs kept: ${kept:-NA}"
                if [ -z "${kept}" ] || (( kept == 0 )); then
                    problem "${t}: no SNPs kept"
                elif [ -n "${prev}" ] && (( kept > prev )); then
                    warn "${t} keeps more SNPs than the looser setting before it"
                fi
                prev=${kept}
            done
        } ;;

    03-select-time/01-baypass-init)
        while IFS=, read -r t k; do
            JOB_IDS+=(${t}-sub${k})
            FILES+=(${S}/dna/baypass/init/baypass-init-${t}-sub${k}.tar.gz)
        done < ${ST}/baypass_init_info.txt ;;

    03-select-time/02-check-omegas)
        JOB_IDS=(check-omegas)
        FILES=(${S}/dna/baypass/init/omega-check.tsv)
        extra () {
            grep -hE "^ *sub[0-9]|OMEGA CHECK" ${LOG_DIR}/check-omegas.out 2> /dev/null
        } ;;

    03-select-time/03-baypass-beta)
        while IFS=, read -r b r seed t; do
            JOB_IDS+=(${t}-beta${b}-run${r})
            FILES+=(${S}/dna/baypass/beta/baypass-${t}-beta${b}-run${r}.tar)
        done < ${ST}/baypass_beta_info.txt ;;

    *)
        echo "Unknown step ${STEP_PATH}" 1>&2
        exit 2 ;;
esac



# =============================================================================
# 1-2. Job status and run time
# =============================================================================

echo "== ${STEP_PATH}: ${#JOB_IDS[@]} expected jobs"

n_ok=0; n_fail=0; n_run=0; n_missing=0
max_h=0; max_id=""
now=$(date +%s)
for id in "${JOB_IDS[@]}"; do
    f=${LOG_DIR}/${id}${LOG_SUFFIX}.log
    if [ ! -f "${f}" ]; then
        n_missing=$((n_missing + 1)); missing_ids+=" ${id}"; continue
    fi
    start=$(grep -E "^001 " "${f}" | tail -n 1)
    end=$(grep -E "^005 " "${f}" | tail -n 1)
    rv=$(grep -oE "return value [0-9]+" "${f}" | tail -n 1 | awk '{print $3}')
    if [ -n "${end}" ] && [ "${rv}" == "0" ]; then
        n_ok=$((n_ok + 1))
    elif [ -n "${end}" ]; then
        n_fail=$((n_fail + 1)); fail_ids+=" ${id}(${rv:-abnormal})"
    else
        n_run=$((n_run + 1)); run_ids+=" ${id}"
    fi
    if [ -n "${start}" ]; then
        t0=$(event_secs "${start}")
        t1=${now}; [ -n "${end}" ] && t1=$(event_secs "${end}")
        if [ -n "${t0}" ] && [ -n "${t1}" ]; then
            h=$(awk -v a=${t0} -v b=${t1} 'BEGIN {printf "%.1f", (b - a) / 3600}')
            if awk -v h=${h} -v m=${max_h} 'BEGIN {exit !(h > m)}'; then
                max_h=${h}; max_id=${id}
            fi
        fi
    fi
done

echo "Finished OK: ${n_ok}; failed: ${n_fail}; running or held: ${n_run}; no log: ${n_missing}"
(( n_fail > 0 )) && problem "failed jobs:${fail_ids}"
(( n_missing > 0 )) && problem "jobs with no log (not submitted?):${missing_ids}"
(( n_run > 0 )) && echo "Running or held:${run_ids}"
if [ -n "${max_id}" ]; then
    echo "Longest run time: ${max_h} h (${max_id})"
    awk -v h=${max_h} -v w=${WARN_H} 'BEGIN {exit !(h > w)}' \
        && warn "a job has run > ${WARN_H} h (CHTC's default limit is 72 h)"
fi



# =============================================================================
# 3. Expected outputs on staging
# =============================================================================

n_files_missing=0
for f in "${FILES[@]}"; do
    if [ ! -f "${f}" ]; then
        n_files_missing=$((n_files_missing + 1))
        (( n_files_missing <= 10 )) && echo "  missing: ${f}"
    fi
done
echo "Outputs on staging: $(( ${#FILES[@]} - n_files_missing )) of ${#FILES[@]}"
if (( n_files_missing > 0 )); then
    if (( n_run > 0 )); then
        warn "${n_files_missing} outputs missing (jobs still running)"
    else
        problem "${n_files_missing} outputs missing"
    fi
fi



# =============================================================================
# 4. Step-specific summary
# =============================================================================

extra

if (( PROBLEMS == 0 )); then
    echo "RESULT: PASS"
    exit 0
else
    echo "RESULT: ${PROBLEMS} problem(s)"
    exit 1
fi
