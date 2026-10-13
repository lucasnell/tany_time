#!/bin/bash

#' Local check (not a cluster job): does lowering npstat's -maxcov to the
#' maximum coverage present change its output, and how much faster is it?
#'
#' Builds npstat from the commit used in the container, simulates a
#' 200-kb pileup (coverage 40-150, ~1% SNPs), and runs npstat with
#' -maxcov 2000 and with -maxcov = observed maximum. Requires gcc and GSL
#' (e.g., `brew install gsl`). Run from the repo root.
#'
#' Output: `_data/npstat-maxcov-check.txt`.
#'

set -e

OUT="$(pwd)/_data/npstat-maxcov-check.txt"
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
    for pos in range(1, 200001):
        if random.random() < 0.05:
            continue
        cov = random.randint(40, 150)
        ref = random.choice(bases)
        if random.random() < 0.01:
            alt = random.choice([b for b in bases if b != ref])
            na = max(1, int(cov * random.uniform(0.05, 0.5)))
            reads = "." * (cov - na) + alt * na
            s.write(f"{pos}\n")
        else:
            reads = "." * cov
        f.write(f"C1\t{pos}\t{ref}\t{cov}\t{reads}\t{'~' * cov}\n")
PYEOF

MX=$(awk 'BEGIN{m=4} {if ($4 > m) m = $4} END{print m}' toy.pileup)
cp toy.pileup big.pileup
cp toy.pileup obs.pileup

t0=$(python3 -c "import time; print(time.time())")
./npstat -n 80 -l 10000 -maxcov 2000 -nolowfreq 0 -snpfile toy.snp \
    big.pileup > /dev/null
t1=$(python3 -c "import time; print(time.time())")
./npstat -n 80 -l 10000 -maxcov ${MX} -nolowfreq 0 -snpfile toy.snp \
    obs.pileup > /dev/null
t2=$(python3 -c "import time; print(time.time())")

if cmp -s big.pileup.stats obs.pileup.stats; then SAME="yes"; else SAME="no"; fi

{
    echo "npstat commit: ${SHA}"
    echo "toy pileup: 200 kb, coverage 40-150, n = 80, 10-kb windows"
    echo "observed max coverage: ${MX}"
    python3 -c "print('seconds, -maxcov 2000: %.2f' % (${t1} - ${t0}))"
    python3 -c "print('seconds, -maxcov ${MX}: %.2f' % (${t2} - ${t1}))"
    echo "identical output: ${SAME}"
} | tee ${OUT}

cd /
rm -rf ${WORK}
