#!/bin/bash

#' Local check (not a cluster job): do the pipeline's tool calls work with
#' the versions pinned in `_bash/_apptainer/*.yml`?
#'
#' Installs the pinned versions (macOS builds; same command-line interface as
#' the container's linux-64 builds) and runs each call with the scripts' own
#' flags on toy data. Where an output could have changed, compares it with
#' the v1.0.0 pins (Python 3.13.5, bcftools 1.22, vcftools 0.1.16).
#' poolfstat is installed from CRAN into a temporary R library because the
#' macOS conda R build segfaults on load.
#' Requires mamba and R. Takes ~15 min. Run from the repo root.
#' Set ENV_DIR to reuse installed environments between runs.
#'
#' Output: `_data/v2-tool-check.txt` (one PASS/FAIL line per check).
#'

ROOT=$(pwd)
OUT=${ROOT}/_data/v2-tool-check.txt
BIN=${ROOT}/_bash/_apptainer/bin
YML=${ROOT}/_bash/_apptainer
WORK=$(mktemp -d)
ENV_DIR=${ENV_DIR:-${WORK}/envs}
mkdir -p ${ENV_DIR}
cd ${WORK}

# Version pinned for a package in an env file, e.g. `pin main-env fastp`
pin () { grep -E "^\s+- $2=" ${YML}/$1.yml | sed 's/.*=//'; }

mk_env () {
    local name=$1; shift
    [ -d ${ENV_DIR}/${name} ] && return 0
    mamba create -y -q -p ${ENV_DIR}/${name} -c conda-forge -c bioconda \
        "$@" > /dev/null 2>&1 \
        || CONDA_SUBDIR=osx-64 mamba create -y -q -p ${ENV_DIR}/${name} \
            -c conda-forge -c bioconda "$@" > /dev/null 2>&1 \
        || echo "could not install $*" 1>&2
}

mk_env py python=$(pin main-env python)
mk_env py_old python=3.13.5
mk_env fastp fastp=$(pin main-env fastp)
mk_env fastqc fastqc=$(pin qc-env fastqc)
mk_env sam samtools=$(pin main-env samtools) bcftools=$(pin main-env bcftools) \
    minimap2=$(pin main-env minimap2) bbmap=$(pin main-env bbmap) \
    vcftools=$(pin main-env vcftools)
mk_env old bcftools=1.22 vcftools=0.1.16
mk_env kr kraken2=$(pin kraken-env kraken2) seqtk=$(pin kraken-env seqtk)
mk_env mito flye=$(pin mito-env flye) medaka=$(pin mito-env medaka) \
    polypolish=$(pin mito-env polypolish)
mk_env bp baypass=$(pin baypass-env baypass)

E=${ENV_DIR}
PY=${E}/py/bin/python
PY0=${E}/py_old/bin/python
P="--remove-indels --max-alleles 2 --window-pi 1000 --window-pi-step 1000"
RES=${WORK}/results.txt
: > ${RES}

# check <description> <command...>: PASS if the command exits 0
check () {
    local d=$1; shift
    if "$@" > /dev/null 2>&1; then echo "PASS: ${d}"; else echo "FAIL: ${d}"; fi \
        | tee -a ${RES}
}
# flags <description> <help text> <flags...>: PASS if all flags are in help
flags () {
    local d=$1 h=$2 o miss=""; shift 2
    for o in "$@"; do grep -q -- "${o}" <<< "${h}" || miss+=" ${o}"; done
    if [ -z "${miss}" ]; then echo "PASS: ${d} flags"
    else echo "FAIL: ${d} flags missing:${miss}"; fi | tee -a ${RES}
}
same () { cmp -s "$1" "$2" && [ -s "$1" ]; }

# ---- toy data
${PY} - << 'PYEOF'
import random, gzip
random.seed(3)
B = "ACGT"
comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
ref = "".join(random.choice(B) for _ in range(5000))
# ref2: SNPs every 97 bp; ref3: ref2 plus SNPs every 211 bp
r2 = list(ref)
for i in range(25, 5000, 97): r2[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[r2[i]]
r3 = list(r2)
for i in range(60, 5000, 211): r3[i] = {"A": "G", "C": "T", "G": "A", "T": "C"}[r3[i]]
for fn, s in [("ref.fa", ref), ("ref2.fa", "".join(r2)), ("ref3.fa", "".join(r3))]:
    open(fn, "w").write(f">CONTIG1\n{s}\n")
# paired 150-bp reads from ref
with gzip.open("r_R1.fastq.gz", "wt") as f1, gzip.open("r_R2.fastq.gz", "wt") as f2:
    for i in range(3000):
        s = random.randint(0, 4600); fr = ref[s:s + random.randint(180, 380)]
        r1 = fr[:150]; r2r = "".join(comp[c] for c in reversed(fr))[:150]
        q = "".join(random.choice("FF:F,") for _ in range(150))
        f1.write(f"@r{i} 1\n{r1}\n+\n{q[:len(r1)]}\n")
        f2.write(f"@r{i} 2\n{r2r}\n+\n{q[:len(r2r)]}\n")
# long reads: half from ref, half from ref2 (heterozygous at ref2's SNPs)
with open("long.fq", "w") as f, open("long_ref.fq", "w") as g:
    for i in range(4000):
        a = random.randint(0, 4000)
        r = (ref if i % 2 else "".join(r2))[a:a + 800]
        f.write(f"@m{i}\n{r}\n+\n{'5' * len(r)}\n")
        r = ref[a:a + 800]
        g.write(f"@n{i}\n{r}\n+\n{'5' * len(r)}\n")
# pileup for the Python SNAPE-input chain
with open("toy.pileup", "w") as p, open("te.gff", "w") as t:
    t.write("CONTIG1\trm\trepeat\t101\t150\t.\t+\t.\t.\n")
    for i, r in enumerate(ref, 1):
        cov = random.randint(10, 80); alt = random.choice([b for b in B if b != r])
        fa = 0.3 if i % 50 == 0 else 0.0
        s = "".join((alt if random.random() < fa else ".") for _ in range(cov))
        if i % 400 == 0: s = s[:-1] + "+2AC" + s[-1]
        p.write(f"CONTIG1\t{i}\t{r}\t{cov}\t{s}\t{'I' * cov}\n")
# 6-pool sync for poolfstat / BayPass
with open("pools.sync", "w") as f:
    for i in range(1, 1501):
        p = random.uniform(0.15, 0.85); row = ["CONTIG1", str(i), "A"]
        for j in range(6):
            n = random.randint(20, 60); q = min(max(p + random.gauss(0, 0.08), 0), 1)
            a = sum(random.random() < q for _ in range(n))
            row.append(f"{a}:{n - a}:0:0:0:0")
        f.write("\t".join(row) + "\n")
PYEOF

# ---- 01-fastp (fastp, FastQC) and pool-qc-files.py
check "fastp $(pin main-env fastp) with 01-fastp flags" \
    ${E}/fastp/bin/fastp --in1 r_R1.fastq.gz --in2 r_R2.fastq.gz \
    --out1 t1.fq.gz --out2 t2.fq.gz --thread 2 --correction \
    --disable_quality_filtering
mkdir -p fqc
check "fastqc $(pin qc-env fastqc) --outdir" \
    ${E}/fastqc/bin/fastqc t1.fq.gz t2.fq.gz --outdir fqc
check "pool-qc-files.py reads fastp JSON (gc and dup not NA)" \
    bash -c "${PY} ${BIN}/pool-qc-files.py --pool x --fastp fastp.json \
        | awk -F'\t' '\$2 != \"NA\" && \$3 != \"NA\" {f = 1} END {exit !f}'"

# ---- 02-bwa / 03-mpileup / 05-kraken samtools calls, bbmerge
export PATH=${E}/sam/bin:${PATH}
samtools faidx ref.fa
gunzip -c t1.fq.gz > u1.fq; gunzip -c t2.fq.gz > u2.fq
check "bbmerge.sh $(pin main-env bbmap)" \
    bbmerge.sh in1=u1.fq in2=u2.fq out=merged.fq outu1=un1.fq outu2=un2.fq
minimap2 -x sr -a ref.fa un1.fq un2.fq 2> /dev/null \
    | samtools view -@ 2 -bh - > all.bam
check "samtools view -f 0x4" samtools view -@ 2 -bh -f 0x4 all.bam -o unm.bam
check "samtools view -q 20 -F 0x100" \
    bash -c "samtools view -@ 2 -bh -q 20 -F 0x100 all.bam > filt.bam"
samtools sort -o sorted.bam filt.bam 2> /dev/null
check "samtools index -@" samtools index -@ 1 sorted.bam
check "samtools idxstats" samtools idxstats sorted.bam
printf "CONTIG1\t0\t5000\n" > t.bed
check "samtools mpileup -B -Q 25 -l -f (non-empty)" \
    bash -c "samtools mpileup -B -Q 25 -l t.bed -f ref.fa sorted.bam \
        2> /dev/null | grep -q ."
check "samtools fastq -@" bash -c "samtools fastq -@ 2 unm.bam > /dev/null"

# ---- 04-ont-pi: minimap2, bcftools, vcftools
for k in 2 3; do
    samtools faidx ref${k}.fa
    minimap2 -x map-ont -a -t 2 -K 1G -2 ref${k}.fa long.fq 2> /dev/null \
        | samtools sort -o long${k}.bam - 2> /dev/null
done
check "minimap2 $(pin main-env minimap2) -x map-ont -K 1G -2 | samtools sort" \
    test -s long3.bam
check "samtools faidx --length 80" samtools faidx --length 80 ref.fa
check "samtools view --threads" \
    bash -c "samtools view -bh -q 20 -F 0x100 --threads 2 long3.bam > /dev/null"
call () {  # call <bcftools dir> <ref> <bam> <-O type> <out>
    $1/bcftools mpileup -Ou --config ont -f $2 $3 2> /dev/null \
        | $1/bcftools call -P 0.01 -mv --threads 2 $4 -o $5 2> /dev/null
}
check "bcftools $(pin main-env bcftools) mpileup --config ont | call -P -mv" \
    call ${E}/sam/bin ref3.fa long3.bam -Oz new.vcf.gz
call ${E}/old/bin ref3.fa long3.bam -Ob old.bcf
${E}/old/bin/vcftools --bcf old.bcf --out old ${P} > /dev/null 2>&1
check "vcftools $(pin main-env vcftools) --gzvcf (04-ont-pi route)" \
    ${E}/sam/bin/vcftools --gzvcf new.vcf.gz --out new ${P}
check "pi identical to v1.0.0 route (bcftools 1.22 BCF, vcftools 0.1.16)" \
    bash -c "cmp -s <(cut -f 1-4,6 new.windowed.pi | tail -n +2) \
        <(tail -n +2 old.windowed.pi) && [ \$(wc -l < old.windowed.pi) -gt 1 ]"
# ---- other envs: flags only
flags "kraken2 $(pin kraken-env kraken2)" "$(${E}/kr/bin/kraken2 --help 2>&1)" \
    --db --threads --confidence --report --output
flags "kraken2-build" "$(${E}/kr/bin/kraken2-build --help 2>&1)" \
    --db --build --threads --clean
flags "seqtk sample" "$(${E}/kr/bin/seqtk sample 2>&1)" "-s"
flags "flye $(pin mito-env flye)" "$(${E}/mito/bin/flye --help 2>&1)" \
    --nano-raw --out-dir --threads --genome-size --asm-coverage
flags "medaka_consensus $(pin mito-env medaka)" \
    "$(PATH=${E}/mito/bin:${PATH} medaka_consensus -h 2>&1)" \
    "-i  " "-d  " "-o  " "-m  " "-t  "
check "medaka model r941_min_hac_g507 exists" bash -c \
    "${E}/mito/bin/medaka tools list_models 2>&1 | grep -q r941_min_hac_g507"
flags "polypolish filter $(pin mito-env polypolish)" \
    "$(${E}/mito/bin/polypolish filter --help 2>&1)" --in1 --in2 --out1 --out2

# ---- Python scripts: identical output under old and new Python
for v in py_old py; do
    mkdir ${v}; (cd ${v}
    X=${E}/${v}/bin/python
    ${X} ${BIN}/PickleRef.py --ref ../ref.fa --output ref
    ${X} ${BIN}/Mpileup2Sync.py --mpileup ../toy.pileup --ref ref.ref \
        --output mp --base-quality-threshold 25 --coding 1.8 --minIndel 5
    ${X} ${BIN}/MaskSYNC_snape.py --sync mp.sync --output masked \
        --indel mp.indel --coverage mp.cov --te ../te.gff --mincov 15 \
        --maxcov 0.95 --maxsnape 0.9
    # add a SNAPE-like 5th column (P(polymorphic) in field 5)
    for x in a b; do
        awk -v s=${x} 'BEGIN {OFS = "\t"; srand(s == "a" ? 1 : 2)}
            $4 !~ /^\./ {pp = ($2 % 50 == 0) ? 0.99 : rand() * 0.5;
                         $5 = "A:T:0.1:0.2:" pp ":0.01:0.0"} 1' \
            masked_masked.sync | gzip -n > ${x}.sync.gz
    done
    ${X} ${BIN}/combine-syncs.py --clean-mono 0.2 --min-pp 0.4 -o counts.sync \
        a.sync.gz b.sync.gz
    ${X} ${BIN}/combine-syncs.py --snape --min-pp 0.4 -o snape.sync \
        a.sync.gz b.sync.gz
    ${X} ${BIN}/pool-qc-sync.py --min-pp 0.4 -o qc.tsv a.sync.gz b.sync.gz
    ) > /dev/null 2>&1
done
PV="Python $(${PY0} -V | cut -d' ' -f2) vs $(${PY} -V | cut -d' ' -f2)"
for f in mp.sync masked_masked.sync masked.bed counts.sync snape.sync qc.tsv; do
    check "${PV}: ${f} identical" same py_old/${f} py/${f}
done
check "${PV}: pickles (ref, cov, indel) load to equal objects" ${PY} -c "
import pickle
for f in ['ref.ref', 'mp.cov', 'mp.indel']:
    assert pickle.load(open('py_old/' + f, 'rb')) == pickle.load(open('py/' + f, 'rb'))"

# ---- 03-select-time: poolfstat and BayPass
PF=$(pin baypass-env r-poolfstat)
mkdir rlib
Rscript -e "
install.packages('https://cloud.r-project.org/src/contrib/Archive/poolfstat/poolfstat_${PF}.tar.gz',
                 lib = 'rlib', repos = NULL, type = 'source', quiet = TRUE)
if (!requireNamespace('poolfstat', lib.loc = 'rlib', quietly = TRUE))
    install.packages('poolfstat', lib = 'rlib', repos = 'https://cloud.r-project.org',
                     quiet = TRUE)" > /dev/null 2>&1
check "poolfstat ${PF}: popsync2pooldata, slots, pooldata2genobaypass" Rscript -e "
suppressPackageStartupMessages(library(poolfstat, lib.loc = 'rlib'))
stopifnot(as.character(packageVersion('poolfstat', lib.loc = 'rlib')) == '${PF}')
pd <- popsync2pooldata(sync.file = 'pools.sync', poolsizes = rep(60, 6),
                       poolnames = paste0('P', 1:6), nthreads = 1)
keep <- rowMeans(pd@readcoverage > 0) >= 0.5
pd@nsnp <- sum(keep)
pd@refallele.readcount <- pd@refallele.readcount[keep, ]
pd@readcoverage <- pd@readcoverage[keep, ]
pd@snp.info <- pd@snp.info[keep, ]
pooldata2genobaypass(pd, prefix = 'toy', subsamplesize = -1L)"
printf "0.1 -0.5 0.3 1.2 -1.0 0.4\n" > toy-log_n.txt
BPV=$(pin baypass-env baypass)
check "BayPass ${BPV} with 01-baypass-init options" \
    ${E}/bp/bin/g_baypass -gfile toy.genobaypass -poolsizefile toy.poolsize \
    -outprefix init -nval 100 -nthreads 2 -seed 115044731
check "BayPass ${BPV} with 03-baypass-beta options" \
    ${E}/bp/bin/g_baypass -gfile toy.genobaypass -poolsizefile toy.poolsize \
    -efile toy-log_n.txt -omegafile init_mat_omega.out -auxmodel \
    -isingbeta 0 -nval 100 -outprefix beta -nthreads 2 -seed 93102805
check "BayPass: files gzipped by 03-baypass-beta exist" test -f \
    beta_summary_pij.out -a -f beta_summary_yij.out \
    -a -f beta_summary_pi_xtx.out -a -f beta_summary_betai.out
check "BayPass: summary_betai.out header as read in _R/" bash -c \
    "[ \"\$(head -1 beta_summary_betai.out | xargs)\" == \
       'COVARIABLE MRK M_Beta SD_Beta PIP BF(dB)' ]"
check "BayPass: summary_pi_xtx.out has M_P, M_XtX, XtXst" bash -c \
    "head -1 beta_summary_pi_xtx.out | tr -s ' ' '\n' \
        | grep -c -x 'M_P\|M_XtX\|XtXst' | grep -q 3"

# ---- write
{
    echo "Tool versions from _bash/_apptainer/*.yml, checked $(date +%Y-%m-%d)"
    echo
    cat ${RES}
    echo
    echo "PASS: $(grep -c '^PASS' ${RES}); FAIL: $(grep -c '^FAIL' ${RES})"
} > ${OUT}
cat ${OUT}

cd /
rm -rf ${WORK}
