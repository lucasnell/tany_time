#!/bin/bash

#' This script calls npstat to estimate population-genetic parameters from
#' masked sync files based on Pool-seq data.
#'
#' Parameters it estimates:
#' - nucleotide diversity
#' - Watterson’s theta
#' - Tajima’s D
#'



. /app/.bashrc
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/lnell/dna/snape"
export PARENT_DIR_OUT="/staging/lnell/dna/npstat"
export GENOME_FULL_PATH="/staging/lnell/Tgraci_assembly.fasta.gz"





#' ========================================================================
#' Inputs
#' ========================================================================

#' First input is read base as usual
export READ_BASE=$1

#' The second input indicates the number of adult midges in this sample.
export N_ADULTS=$2

#' The third indicates the size of sliding window in kb.
export WIN_KB=$3

export IN_SNAPE_PREFIX=${READ_BASE}_snape
export SNP=${IN_SNAPE_PREFIX}_masked.snp
export SYNC=${IN_SNAPE_PREFIX}_part_masked.sync.gz

export GENOME=$(basename ${GENOME_FULL_PATH%.gz})


if [ ! -f ${PARENT_DIR_IN}/${IN_SNAPE_PREFIX}.tar ]; then
    echo "${PARENT_DIR_IN}/${IN_SNAPE_PREFIX}.tar does not exist! " 1>&2
    safe_exit 111
fi
if [ ! -f ${GENOME_FULL_PATH} ]; then
    echo "${GENOME_FULL_PATH} does not exist! " 1>&2
    safe_exit 222
fi
# Contig names in the same order as the reference.
# Used to organize output files since SNAPE analyzes by contig.
# Make sure contig names don't have spaces in them!
export CONTIG_NAMES=($(gunzip -c ${GENOME_FULL_PATH} | grep "^>" | sed 's/>//g' | sed 's/\s.*$//'))
check_exit_status "get contig names" $?




#' ========================================================================
#' Outputs
#' ========================================================================

# Final files / directories
export OUT_DIR=${READ_BASE}_npstat_${WIN_KB}kb
export OUT_FILE=${OUT_DIR}.stat
# Intermediate file:
export PILEUP=${READ_BASE}_snape_part_masked.pileup




#' ========================================================================
#' Copy over inputs
#' ========================================================================

mkdir ${OUT_DIR}

cd ${OUT_DIR}

tar -xf ${PARENT_DIR_IN}/${IN_SNAPE_PREFIX}.tar ${IN_SNAPE_PREFIX}/${SNP}.gz && \
    gunzip ${IN_SNAPE_PREFIX}/${SNP}.gz && \
    mv ${IN_SNAPE_PREFIX}/${SNP} ./
check_exit_status "extract SNP" $?

tar -xf ${PARENT_DIR_IN}/${IN_SNAPE_PREFIX}.tar ${IN_SNAPE_PREFIX}/${SYNC} && \
    mv ${IN_SNAPE_PREFIX}/${SYNC} ./
check_exit_status "extract SYNC" $?

rm -r ${IN_SNAPE_PREFIX}



#' ========================================================================
#' Convert sync file to pileup
#' ========================================================================

# Takes ~5 min

python3 << EOF
"""
Convert gSYNC files (optionally gzipped) into a pileup file.

NOTE: does not work with indels!

"""

import sys
import os.path
import gzip
from datetime import datetime

if __name__ == "__main__":

    out = "${PILEUP}"
    sync_in = "${SYNC}"

    if not sync_in.endswith((".sync", ".sync.gz")):
        print("Strange suffix to input file name (must be .sync or " +
              ".sync.gz). Exiting.", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(sync_in):
        print(sync_in + " does not exist. Exiting.", file=sys.stderr)
        sys.exit(1)
    if not out.endswith((".pileup", ".pileup.gz", ".txt", ".txt.gz")):
        print("Strange suffix to output file name (must be .pileup, " +
              ".pileup.gz, .txt, or .txt.gz). Exiting.", file=sys.stderr)
        sys.exit(1)
    out_dir_file = os.path.split(out)
    if out_dir_file[0] != "" and not os.path.exists(out_dir_file[0]):
        print("output dir (" + out_dir_file[0] + ") not found", file=sys.stderr)
        sys.exit(1)
    if os.path.exists(out):
        print(out + " already exists. Overwriting not allowed.\n", file=sys.stderr)
        sys.exit(1)

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print("Converting sync file to " + out)

    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("started reading/writing (" + time_str + ")")

    if sync_in.endswith(".gz"):
        in_file = gzip.open(sync_in, "rt")
    else:
        in_file = open(sync_in, "r")

    if out.endswith(".gz"):
        out_file = gzip.open(out, "wt")
    else:
        out_file = open(out, "w")

    # Bases in the order they're listed in the sync file's 4th column.
    bases = "ATCGN"
    # Vector to map characters to 0-based indices for where the allele
    # count should be in the sync file's 4th column.
    base_map = [-1 for x in range(255)]
    for i, x in enumerate(bases):
        base_map[ord(x)] = i
        base_map[ord(x.lower())] = i


    # this gets changed each iteration to have "." used for reference:
    tmp_bases = [x for x in bases]

    for i, line in enumerate(in_file):
        l_split = line.rstrip().split("\t")
        if len(l_split) < 4:
            print("only " + str(len(l_split)) + " columns in file " +
                  sync_in + " line " + str(i), file=sys.stderr)
            sys.exit(1)
        if l_split[3] == ".:.:.:.:.:.":
            continue
        new_line = "\t".join(l_split[:3])
        ref = l_split[2]
        ref_idx = base_map[ord(ref)]
        tmp_bases[ref_idx] = "."
        counts = [int(x) for x in l_split[3].split(":")]
        new_line += "\t"
        new_line += str(sum(counts))
        new_line += "\t"
        for j in range(len(tmp_bases)):
            new_line += tmp_bases[j] * counts[j]
        new_line += "\t"
        new_line += "~" * sum(counts)
        new_line += "\n"
        b = out_file.write(new_line)
        j += 1
        tmp_bases[ref_idx] = bases[ref_idx]

    in_file.close()
    out_file.close()

    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("... finished (" + time_str + ")")

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    sys.exit(0)

EOF


rm ${SYNC}




#' ========================================================================
#' Split pileup and SNP files by contig
#' ========================================================================
mkdir tmp
cd tmp

# Pileup file first:
# This does the splitting:
# Takes ~1 min on a 23G pileup
awk '{if (last != $1) close(last); print >> $1; last = $1}' ../${PILEUP}
# Rename contig files and verify they exist in the reference
for contig in *; do
    if [[ ! " ${CONTIG_NAMES[*]} " =~ " ${contig} " ]]; then
        echo "Contig ${contig} does not exist in reference! " 1>&2
        if [[ $- != *i* ]]; then
            exit 1
        else
            break
        fi
    fi
    mv $contig ../${contig}.pileup
done
rm ../${PILEUP}


# Now same for SNP file
awk '{if (last != $1) close(last); print $2 >> $1; last = $1}' ../${SNP}
for contig in *; do
    if [[ ! " ${CONTIG_NAMES[*]} " =~ " ${contig} " ]]; then
        echo "Contig ${contig} does not exist in reference! " 1>&2
        if [[ $- != *i* ]]; then
            exit 1
        else
            break
        fi
    fi
    mv $contig ../${contig}.snp
done
rm ../${SNP}



# Now go back to main directory
cd ..
rm -r tmp

# # If you want to check which files are present
# for contig in ${CONTIG_NAMES[@]}; do
#     if [ ! -f ${contig}.pileup ]; then
#         echo "${contig}.pileup not found"
#     fi
#     if [ ! -f ${contig}.snp ]; then
#         echo "${contig}.snp not found"
#     fi
# done


#' ========================================================================
#' npstat
#' ========================================================================

# Shell script to run npstat on a single contig.
# Each one takes a while to run, so I want to do this in parallel.
cat << EOF > one_npstat.sh
#!/bin/bash
export contig=\${1}
if [ ! -f \${contig}.pileup ] || [ ! -f \${contig}.snp ]; then
    rm \${contig}.pileup \${contig}.snp 2> /dev/null
    exit 1
fi
npstat -n $((N_ADULTS * 2)) -l ${WIN_KB}000 -maxcov 10000 -nolowfreq 0 \\
    -snpfile \${contig}.snp \${contig}.pileup \\
    &> \${contig}_npstat.out
status=\$?
if [ "\$status" != "0" ]; then
    exit 2
fi
rm \${contig}.pileup \${contig}.snp \${contig}_npstat.out

# Keep only columns that are useful given my input parameters
cut -f 1,2,4-10 \${contig}.pileup.stats \\
    > \${contig}.stats
rm \${contig}.pileup.stats

exit 0
EOF

chmod +x one_npstat.sh


cat << EOF > all_npstats.py
#!/usr/bin/env python3
import subprocess as sp
import sys

def work(contig):
    """Defines the work unit on an input file"""
    cmd = "./one_npstat.sh " + contig
    ret = sp.run(cmd, shell = True)
    return ret

if __name__ == "__main__":
    #
    from dask.distributed import Client, LocalCluster
    from dask import config as dc
    dc.set({"distributed.admin.system-monitor.gil.enabled": False})
    ##> # for interactive job where you want progress bar:
    ##> from dask.distributed import progress
    #
    tasks = "${CONTIG_NAMES[@]}".split(" ")
    n_tasks = len(tasks)
    #
    print("Starting npstat...")
    with LocalCluster(processes=True,
        n_workers=${THREADS},
        threads_per_worker=1) as cluster, Client(cluster) as client:
        futures = []
        for task in tasks:
            future = client.submit(work, task)
            futures.append(future)
        ##> # for interactive job where you want progress bar:
        ##> progress(futures)
        #
        # Wait until futures done, then collect output:
        results = client.gather(futures)
    #
    print("... Ended npstat")
    # Now just get the return codes:
    return_codes = [x.returncode for x in results]
    #
    for i in range(n_tasks):
        if return_codes[i] == 1:
            print("contig " + tasks[i] + " not found!", file=sys.stderr)
        elif return_codes[i] == 2:
            print("npstat failed on contig " + tasks[i] + "!", file=sys.stderr)
            print("See " + tasks[i] + "_npstat.out", file=sys.stderr)
    sys.exit(0)
    #
EOF

chmod +x all_npstats.py

# If it's an interactive job, uncomment lines for progress bar:
if [[ $- == *i* ]]; then sed -i 's/##> //g' all_npstats.py; fi

# Takes ~1 hr with 16 threads
./all_npstats.py
status=$?

rm all_npstats.py one_npstat.sh



# Combine output files.
# Add header from one output file (all are the same) to the final output.
echo -ne "sequence\t" > ${OUT_FILE}
head -n 1 ${CONTIG_NAMES[0]}.stats >> ${OUT_FILE}
# I'm using `CONTIG_NAMES` so that the combined file is in the same order
# as the reference.
for contig in ${CONTIG_NAMES[@]}; do
    # remove header, add the contig column, then add this to main output file
    tail -n+2 ${contig}.stats \
        | awk -v CONTIG_NAME="${contig}" '{print CONTIG_NAME"\t"$0}' - \
        >> ${OUT_FILE} && \
        rm ${contig}.stats
done




#' Modified by LAN on 2022-04-28 bc output from `npstat` had columns 9-10 and
#' 11-12 switched compared to `npstat` manual
# 1. window number,
# 2. number of bases covered by sequences,
# 3. number of bases covered and with known outgroup allele,
# 4. average read depth,
# 5. number of segregating sites S,
# 6. Watterson estimator of θ,
# 7. Tajima’s Π estimator of heterozygosity,
# 8. Tajima’s D,
# 9. variance of the number of segregating sites,
# 10. variance of the Watterson estimator,
# 11. unnormalized Fay and Wu’s H,
# 12. normalized Fay and Wu’s H,
# 13. divergence per base (from outgroup),
# 14. nonsynonimous polymorphisms,
# 15. synonimous polymorphisms,
# 16. nonsynonimous divergence,
# 17. synonimous divergence,
# 18. α (fraction of substitutions fixed by positive selection).





#' ========================================================================
#' Handle output file
#' ========================================================================


gzip ${OUT_FILE}


mv ${OUT_FILE}.gz ../
cd ..
# Copying, not moving so that they will be moved to ResearchDrive after exit
cp ${OUT_FILE}.gz ${PARENT_DIR_OUT}/

rm -r ${OUT_DIR}

#' Do NOT use `safe_exit` because that'll remove the files before they
#' can be moved to ResearchDrive.
exit 0
