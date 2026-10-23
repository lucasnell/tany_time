#!/bin/bash

#' Local check (not a cluster job) of the combined Dorado 0.9.6 sup
#' basecalls from `00-basecall/03-combine`, copied from staging to the
#' external drive. Checks that the gzipped FASTQ decompresses fully, holds
#' exactly the summary's reads with mean q >= 10 (no duplicate IDs, sequence
#' and quality lengths equal), and compares read counts, yield, read N50,
#' and predicted quality with the Guppy 5.0.11 hac basecalls, matching reads
#' by ID. Qualities are the basecallers' own predictions (mean_qscore), not
#' measured accuracy. Takes ~30-60 min (reads ~70 GB of gzip). Run from the
#' repo root.
#' The Guppy 5.0.11 hac files it reads were deleted from the local disk
#' on 2026-10-23 (copies kept elsewhere by the owner); restore them to
#' rerun.
#'
#' Output: `_data/dorado-sup-check.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/dorado-sup-check.txt
D_FQ=/Volumes/lan_images4/ont/basecalls_dorado-0.9.6_sup.fastq.gz
D_SUMM=/Volumes/lan_images4/ont/basecalls_dorado-0.9.6_sup--summary.tsv.gz
G_FQ=${HOME}/_data/midgenomes/_sra/basecalls_guppy-5.0.11.fastq.gz
G_SUMM=${HOME}/_data/ont/basecalls_guppy-5.0.11--sequencing_summary.txt.gz
MIN_QS=10
WORK=$(mktemp -d)
cd ${WORK}
export LC_ALL=C

for f in ${D_FQ} ${D_SUMM} ${G_FQ} ${G_SUMM}; do
    [ -s ${f} ] || { echo "missing: ${f}" 1>&2; exit 1; }
done

# read_id, length, mean q from a summary (columns found by name)
summ_cols() {
    gunzip -c $1 | awk -F'\t' -v OFS='\t' '
        NR == 1 {for (i = 1; i <= NF; i++) {if ($i == "read_id") I = i
                     if ($i == "sequence_length_template") L = i
                     if ($i == "mean_qscore_template") Q = i}; next}
        {print $I, $L, $Q}' | sort -k 1,1
}
# lengths on stdin -> "reads bp N50"
len_stats() {
    sort -rn | awk '{l[NR] = $1; t += $1}
        END {for (i = 1; i <= NR; i++) {c += l[i]; if (c >= t / 2) {n50 = l[i]; break}}
             printf "%d reads, %.3f Gb, N50 %d", NR, t / 1e9, n50}'
}
# FASTQ on stdin -> id, sequence length, qs tag (NA if absent); flags
# records whose sequence and quality lengths differ or are malformed
fq_cols() {
    awk -v OFS='\t' '
        NR % 4 == 1 {if (substr($1, 1, 1) != "@") bad++; id = substr($1, 2); q = "NA"
                     for (i = 2; i <= NF; i++) if ($i ~ /^qs:f:/) q = substr($i, 6)}
        NR % 4 == 2 {sl = length($0)}
        NR % 4 == 3 {if (substr($0, 1, 1) != "+") bad++}
        NR % 4 == 0 {if (length($0) != sl) bad++; print id, sl, q}
        END {if (NR % 4 != 0) bad++; print bad + 0 > "/dev/stderr"}'
}

summ_cols ${D_SUMM} > d_summ.tsv
summ_cols ${G_SUMM} > g_summ.tsv

gunzip -c ${D_FQ} | fq_cols 2> d_fq.bad > d_fq_unsorted.tsv
D_GZ=${PIPESTATUS[0]}
sort -k 1,1 d_fq_unsorted.tsv > d_fq.tsv
gunzip -c ${G_FQ} | fq_cols 2> g_fq.bad > g_fq_unsorted.tsv
G_GZ=${PIPESTATUS[0]}
sort -k 1,1 g_fq_unsorted.tsv > g_fq.tsv

{
echo "Dorado 0.9.6 sup (00-basecall/03-combine) vs Guppy 5.0.11 hac; checked $(date +%F)"
echo "files: ${D_FQ}, ${D_SUMM}; ${G_FQ}, ${G_SUMM}"
echo
echo "== Dorado sup FASTQ integrity =="
echo "gunzip exit status: ${D_GZ} (0 = whole file decompressed)"
echo "malformed records or sequence/quality length mismatches: $(cat d_fq.bad)"
echo "duplicate read IDs: $(cut -f 1 d_fq.tsv | uniq -d | wc -l | tr -d ' ')"
echo "reads without qs tag: $(awk '$3 == "NA"' d_fq.tsv | wc -l | tr -d ' ')"
echo "reads with qs < ${MIN_QS}: $(awk -v m=${MIN_QS} '$3 != "NA" && $3 + 0 < m' d_fq.tsv | wc -l | tr -d ' ')"
awk -v m=${MIN_QS} '$3 >= m' d_summ.tsv | cut -f 1 > d_summ_pass.ids
cut -f 1 d_fq.tsv > d_fq.ids
echo "summary reads with q >= ${MIN_QS} missing from FASTQ: $(comm -23 d_summ_pass.ids d_fq.ids | wc -l | tr -d ' ')"
echo "FASTQ reads not in summary at q >= ${MIN_QS}: $(comm -13 d_summ_pass.ids d_fq.ids | wc -l | tr -d ' ')"
join -t $'\t' d_fq.tsv d_summ.tsv \
    | awk '$2 != $4 {n++} END {print "FASTQ length != summary length: " n + 0}'
echo
echo "== Guppy hac FASTQ (used by the mitogenome tests and 04-ont-pi) =="
echo "gunzip exit status: ${G_GZ}; malformed records: $(cat g_fq.bad)"
echo "duplicate read IDs: $(cut -f 1 g_fq.tsv | uniq -d | wc -l | tr -d ' ')"
cut -f 1 g_fq.tsv > g_fq.ids
echo "FASTQ reads not in Guppy summary: $(cut -f 1 g_summ.tsv | comm -13 - g_fq.ids | wc -l | tr -d ' ')"
join -t $'\t' g_fq.tsv g_summ.tsv | cut -f 5 | sort -g \
    | awk '{q[NR] = $1} END {printf "Guppy summary q of FASTQ reads: min %.2f, median %.2f\n", q[1], q[int((NR + 1) / 2)]}'
echo
echo "== Yield =="
printf "Dorado summary, all reads:     %s\n" "$(cut -f 2 d_summ.tsv | len_stats)"
printf "Dorado FASTQ (q >= %s):        %s\n" ${MIN_QS} "$(cut -f 2 d_fq.tsv | len_stats)"
printf "Guppy summary, all reads:      %s\n" "$(cut -f 2 g_summ.tsv | len_stats)"
printf "Guppy summary, q >= %s:        %s\n" ${MIN_QS} "$(awk -v m=${MIN_QS} '$3 >= m' g_summ.tsv | cut -f 2 | len_stats)"
printf "Guppy FASTQ:                   %s\n" "$(cut -f 2 g_fq.tsv | len_stats)"
echo
echo "== Read sets =="
cut -f 1 d_summ.tsv > d_summ.ids
cut -f 1 g_summ.tsv > g_summ.ids
echo "duplicate IDs in Dorado summary: $(uniq -d d_summ.ids | wc -l | tr -d ' ')"
echo "reads in both summaries: $(comm -12 d_summ.ids g_summ.ids | wc -l | tr -d ' ')"
echo "Dorado only: $(comm -23 d_summ.ids g_summ.ids | wc -l | tr -d ' '); Guppy only: $(comm -13 d_summ.ids g_summ.ids | wc -l | tr -d ' ')"
echo "Guppy FASTQ reads in Dorado FASTQ: $(comm -12 g_fq.ids d_fq.ids | wc -l | tr -d ' ') of $(wc -l < g_fq.ids | tr -d ' ')"
echo "Dorado FASTQ reads not in Guppy FASTQ: $(comm -13 g_fq.ids d_fq.ids | wc -l | tr -d ' ')"
echo
echo "== Predicted quality, reads in both summaries =="
join -t $'\t' d_summ.tsv g_summ.tsv > both.tsv
awk '{d += $3; g += $5; if ($2 != $4) nl++}
    END {printf "mean q: Dorado %.2f, Guppy %.2f; length differs in %d of %d reads\n", d / NR, g / NR, nl, NR}' both.tsv
awk '{print $3 - $5}' both.tsv | sort -g | awk '{x[NR] = $1}
    END {printf "q difference (Dorado - Guppy): 10th percentile %.2f, median %.2f, 90th %.2f\n",
         x[int(NR * 0.1) + 1], x[int((NR + 1) / 2)], x[int(NR * 0.9)]}'
for b in "0 1000" "1000 5000" "5000 20000" "20000 1e12"; do
    set -- ${b}
    awk -v lo=$1 -v hi=$2 '$4 >= lo && $4 < hi {n++; d += $3; g += $5}
        END {if (n) printf "Guppy length [%s, %s): %d reads, mean q Dorado %.2f, Guppy %.2f\n", lo, hi, n, d / n, g / n}' both.tsv
done
echo "q bins (reads, all): q < 7 / 7-10 / 10-15 / 15-20 / >= 20"
for f in d_summ.tsv g_summ.tsv; do
    awk -v f=${f%%_*} '{if ($3 < 7) a++; else if ($3 < 10) b++; else if ($3 < 15) c++; else if ($3 < 20) e++; else h++}
        END {printf "  %s: %d / %d / %d / %d / %d\n", (f == "d" ? "Dorado" : "Guppy"), a, b, c, e, h}' ${f}
done
} > ${OUT}

sed -i '' 's/[[:space:]]*$//' ${OUT}
echo "work dir (delete when done): ${WORK}"
