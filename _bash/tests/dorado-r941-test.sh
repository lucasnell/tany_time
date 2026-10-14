#!/bin/bash

#' Local check (not a cluster job): does Dorado 0.9.6 basecall the 2019
#' R9.4.1 FAST5 data with the sup model, and how do its predicted read
#' qualities compare with the Guppy 5.0.11 hac basecalls of the same reads?
#'
#' Downloads Dorado 0.9.6 (Apple Silicon build) and `dna_r9.4.1_e8_sup@v3.6`,
#' extracts the first FAST5 file from the raw-signal archive on the external
#' drive, basecalls N_READS reads, and joins them to the Guppy sequencing
#' summary by read ID. Qualities are the basecallers' own predictions
#' (mean_qscore), not measured accuracy. Requires internet and the external
#' drive. Takes ~10 min. Run from the repo root.
#'
#' Output: `_data/dorado-r941-test.txt`.
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/dorado-r941-test.txt
ARCHIVE=/Volumes/TanyGen5TB/nanopore/ives_fast5.tar.bz2
SUMM=${HOME}/_data/ont/basecalls_guppy-5.0.11--sequencing_summary.txt.gz
DORADO=dorado-0.9.6-osx-arm64
MODEL=dna_r9.4.1_e8_sup@v3.6
N_READS=200
WORK=$(mktemp -d)
cd ${WORK}

[ -f ${ARCHIVE} ] || { echo "external drive not mounted" 1>&2; exit 1; }

curl -sSL -o d.zip https://cdn.oxfordnanoportal.com/software/analysis/${DORADO}.zip \
    && unzip -q d.zip || { echo "Dorado download failed" 1>&2; exit 1; }
D=${WORK}/${DORADO}/bin/dorado
${D} download --model ${MODEL} > /dev/null 2>&1 \
    || { echo "model download failed" 1>&2; exit 1; }

mkdir f5
python3 -I -u - ${ARCHIVE} << 'PYEOF'
import bz2, shutil, sys, tarfile
t = tarfile.open(fileobj=bz2.open(sys.argv[1], "rb"), mode="r|")
for m in t:
    if m.isfile() and m.name.endswith(".fast5"):
        with open("f5/first.fast5", "wb") as o:
            shutil.copyfileobj(t.extractfile(m), o)
        print(m.name)
        break
PYEOF

${D} basecaller --max-reads ${N_READS} --emit-fastq ${MODEL} f5/ \
    > sup.fastq 2> bc.log || { tail bc.log 1>&2; exit 1; }

awk 'NR % 4 == 1 {id = substr($1, 2)
                  for (i = 2; i <= NF; i++) if ($i ~ /^qs:f:/) q = substr($i, 6)}
     NR % 4 == 2 {print id, length($0), q}' sup.fastq | sort -k1,1 > sup.txt
gunzip -c ${SUMM} | awk -F'\t' 'NR > 1 {print $2, $14, $15}' | sort -k1,1 \
    | join sup.txt - > both.txt

{
    echo "Dorado $(${D} --version 2>&1 | tail -n 1), model ${MODEL};" \
        "checked $(date +%Y-%m-%d)"
    echo "first FAST5 file of the archive; ${N_READS} reads requested;" \
        "$(grep -h "Simplex reads basecalled" bc.log | grep -oE "[0-9]+$") basecalled;" \
        "$(wc -l < both.txt | tr -d ' ') matched to the Guppy summary"
    awk '{ds += $3; gs += $5; dl += $2; gl += $4; n++}
         END {printf "mean predicted q: Dorado sup %.2f, Guppy hac %.2f\n", ds / n, gs / n
              printf "mean read length: Dorado %.0f, Guppy %.0f\n", dl / n, gl / n}' both.txt
    awk '{print $3 - $5}' both.txt | sort -n \
        | awk '{a[NR] = $1} END {printf "q difference (sup - hac): median %.2f," \
               " 10th percentile %.2f, 90th %.2f\n", a[int(NR / 2)],
               a[int(NR * 0.1)], a[int(NR * 0.9)]}'
    grep -h "Samples/s" bc.log | sed 's/^.*> /Mac speed: /'
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
