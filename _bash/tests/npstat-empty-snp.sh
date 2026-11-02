#!/bin/bash

#' Local check (not a cluster job): what does npstat do with an empty
#' `-snpfile`? In `02-allele-freqs/02-npstat`, a contig with no SNPs in a
#' pool has no SNP file and gets no windows, whereas a window with no SNPs
#' on a contig with SNPs is reported with S = 0.
#'
#' Builds npstat from the commit used in the container and simulates a
#' 50-kb pileup (coverage 40-150) with variant reads at ~1% of sites,
#' except in the third 10-kb window. Runs npstat with a SNP file listing
#' the variant sites and with an empty SNP file, and compares the windows.
#' Requires gcc and GSL (e.g., `brew install gsl`). Run from the repo root.
#'
#' Output: `_data/npstat-empty-snp-check.txt`.
#'

set -e

OUT="$(pwd)/_data/npstat-empty-snp-check.txt"
SHA="dc230fba3ac434ed6526ab7cd63a67d69d34b7d6"
WORK=$(mktemp -d)
cd ${WORK}

curl -sL -o npstat.zip https://github.com/lucaferretti/npstat/archive/${SHA}.zip
unzip -q npstat.zip
gcc -O2 -o npstat npstat-${SHA}/NPStat-v1.c $(gsl-config --cflags) \
    $(gsl-config --libs) -lm 2> /dev/null

python3 - << 'PYEOF'
import random
random.seed(1)
bases = "ACGT"
with open("toy.pileup", "w") as f, open("toy.snp", "w") as s:
    for pos in range(1, 50001):
        if random.random() < 0.05:
            continue
        cov = random.randint(40, 150)
        ref = random.choice(bases)
        if random.random() < 0.01 and not 20000 < pos <= 30000:
            alt = random.choice([b for b in bases if b != ref])
            na = max(1, int(cov * random.uniform(0.05, 0.5)))
            reads = "." * (cov - na) + alt * na
            s.write(f"{pos}\n")
        else:
            reads = "." * cov
        f.write(f"C1\t{pos}\t{ref}\t{cov}\t{reads}\t{'~' * cov}\n")
PYEOF

: > empty.snp
cp toy.pileup snps.pileup
cp toy.pileup empty.pileup
MX=$(awk 'BEGIN{m=4} {if ($4 > m) m = $4} END{print m}' toy.pileup)

# Same call as `02-npstat.sh` (n = 2 x adults; 80 here):
set +e
./npstat -n 80 -l 10000 -maxcov ${MX} -nolowfreq 0 -snpfile toy.snp \
    snps.pileup > snps.log 2>&1
RC_SNPS=$?
./npstat -n 80 -l 10000 -maxcov ${MX} -nolowfreq 0 -snpfile empty.snp \
    empty.pileup > empty.log 2>&1
RC_EMPTY=$?
set -e

{
    echo "npstat commit: ${SHA}"
    echo "toy pileup: 50 kb, coverage 40-150, n = 80, 10-kb windows;"
    echo "$(wc -l < toy.snp | tr -d ' ') variant sites, none in window 3"
    echo "exit status, SNP file: ${RC_SNPS}; empty SNP file: ${RC_EMPTY}"
    for k in snps empty; do
        echo
        echo "== ${k} (columns kept by 02-npstat.sh)"
        if [ -f ${k}.pileup.stats ]; then
            cut -f 1,2,4-10 ${k}.pileup.stats
        else
            echo "no output; npstat said:"
            cat ${k}.log
        fi
    done
    if [ -f snps.pileup.stats ] && [ -f empty.pileup.stats ]; then
        echo
        if cmp -s <(sed -n 4p snps.pileup.stats) <(sed -n 4p empty.pileup.stats); then
            echo "window 3 (no SNPs either way) identical in both runs: yes"
        else
            echo "window 3 (no SNPs either way) identical in both runs: no"
        fi
        if cmp -s <(cut -f 1-4 snps.pileup.stats) <(cut -f 1-4 empty.pileup.stats); then
            echo "window, length, outgroup length, and read depth identical in all windows: yes"
        else
            echo "window, length, outgroup length, and read depth identical in all windows: no"
        fi
    fi
} | tee ${OUT}

echo "(scratch folder left at ${WORK})"
