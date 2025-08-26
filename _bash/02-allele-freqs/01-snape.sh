#!/bin/bash

#' This script calls snape-pooled to estimate allele frequencies for
#' Pool-seq data.

#' Note #1:
#' For many of the scripts inside `snape-files`, the `output`
#' argument is supposed to be a prefix only.
#' Hence, there are multiple instances below where I remove the extension
#' to the final file names for input to this argument.

#' Note #2:
#' To compress the `snape-files` folder into a tar file that doesn't generate
#' warnings on linux, run the following:
#'
#' tar -cz --no-xattrs --exclude ".*" -f snape-files.tar.gz snape-files
#'


. /app/.bashrc
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/lnell/dna/mpileup"
export PARENT_DIR_OUT="/staging/lnell/dna/snape"
export GENOME_FULL_PATH="/staging/lnell/Tgraci_assembly.fasta.gz"
export REPEATS_FULL_PATH="/staging/lnell/Tgraci_repeats_locs.gff3.gz"


# export READ_BASE=Ash-19_S5
# export N_ADULTS=30



#' ========================================================================
#' Inputs
#' ========================================================================

#' First input is read base as usual
export READ_BASE=$1
#' The second input indicates the number of adult midges in this sample.
export N_ADULTS=$2

export IN_PILEUP=${READ_BASE}_mpileup.txt
export GENOME=$(basename ${GENOME_FULL_PATH%.gz})
export REPEATS=$(basename ${REPEATS_FULL_PATH%.gz})

if [ ! -f ${PARENT_DIR_IN}/${IN_PILEUP}.gz ]; then
    echo "${PARENT_DIR_IN}/${IN_PILEUP} does not exist! " 1>&2
    # Don't actually exit if it's an interactive job:
    if [[ $- != *i* ]]; then exit 111; fi
fi
if [ ! -f ${GENOME_FULL_PATH} ]; then
    echo "${GENOME_FULL_PATH} does not exist! " 1>&2
    if [[ $- != *i* ]]; then exit 222; fi
fi
if [ ! -f ${REPEATS_FULL_PATH} ]; then
    echo "${REPEATS_FULL_PATH} does not exist! " 1>&2
    if [[ $- != *i* ]]; then exit 222; fi
fi



#' ========================================================================
#' Output names and parameter values
#' ========================================================================


# Final files / directories
export OUT_DIR=${READ_BASE}_snape
export UNMASKED_PREFIX=${READ_BASE}_snape
export MASKED_FULL_PREFIX=${READ_BASE}_snape_masked
export MASKED_PART_PREFIX=${READ_BASE}_snape_part_masked
# Intermediates
export MP_INFO_PREFIX=${IN_PILEUP%.txt}_info
export GENOME_PICKLE=${GENOME%.fasta}_pickle.ref


# Set parameters for snape:
export THETA=0.0088
export D=0.01
export PRIOR_TYPE="informative"
export FOLD="unfolded"

# Set parameters for processing snape output:
export MAX_SNAPE=0.9
export MAX_COV=0.95
export MIN_COV=10




#' ========================================================================
#' Prep for downstream steps.
#' ========================================================================

mkdir ${OUT_DIR}
cd ${OUT_DIR}

gunzip -c ${PARENT_DIR_IN}/${IN_PILEUP}.gz \
    > ${IN_PILEUP}
check_exit_status "gunzip ${IN_PILEUP}.gz" $?

# For progress bar in interactive jobs:
# pv ${PARENT_DIR_IN}/${IN_PILEUP}.gz \
#     | gunzip \
#     > ${IN_PILEUP}
# check_exit_status "gunzip ${IN_PILEUP}.gz" $?

gunzip -c ${GENOME_FULL_PATH} \
    > ${GENOME}
check_exit_status "gunzip ${GENOME}.gz" $?

gunzip -c ${REPEATS_FULL_PATH} \
    > ${REPEATS}
check_exit_status "gunzip ${REPEATS}.gz" $?

# Files to process SNAPE output:
tar -xzf ${PARENT_DIR_OUT}/snape-files.tar.gz -C ./
check_exit_status "cp, extract snape-files.tar.gz" $?

# "Pickle" the reference
python3 ./snape-files/PickleRef.py --ref ${GENOME} \
    --output ${GENOME_PICKLE%.ref}
check_exit_status "PickleRef" $?

# Scaffold names in the same order as the reference.
# Used to organize output files since SNAPE analyzes by contig.
# Make sure contig names don't have spaces in them!
export CONTIG_NAMES=($(grep "^>" ${GENOME} | sed 's/>//g' | sed 's/\s.*$//'))



#' ========================================================================
#' Summarize some mpileup file info for use later
#' Produces the following files:
#' - ${MP_INFO_PREFIX}.sync
#' - ${MP_INFO_PREFIX}.cov
#' - ${MP_INFO_PREFIX}.indel
#'
#'
#' Takes ~110 min
python3 ./snape-files/Mpileup2Sync.py \
  --mpileup ${IN_PILEUP} \
  --ref ${GENOME_PICKLE} \
  --output ${MP_INFO_PREFIX} \
  --base-quality-threshold 25 \
  --coding 1.8 \
  --minIndel 5
check_exit_status "Mpileup2Sync" $?

# Not needed hereafter:
rm ${MP_INFO_PREFIX}.sync






#' ========================================================================
#' Split mpileup file by contig
mkdir tmp
cd tmp

# This does the splitting (takes ~14 min):
awk '{if (last != $1) close(last); print >> $1; last = $1}' ../${IN_PILEUP}

# Rename contig files and verify they exist in the reference
for contig in *; do
    if [[ ! " ${CONTIG_NAMES[*]} " =~ " ${contig} " ]]; then
        echo "Scaffold ${contig} does not exist in reference! " 1>&2
        if [[ $- != *i* ]]; then
            exit 1
        else
            break
        fi
    fi
    mv $contig ${contig}_mp.txt
done


# Now go back to main directory
mv *_mp.txt ../
cd ..
rm -r tmp ${IN_PILEUP}




#' ========================================================================
#' SNAPE-pooled


#' SNAPE-pooled can't use multithreading itself, so using python to do
#' it for us.

# Script to run snap-pooled on one contig:
cat << EOF > one_snape.sh
#!/bin/bash
contig=\$1

# If this contig isn't present, we add a filler line to the snape
# output file to make sure that all of our eventual gSYNC files have the
# same number of rows.
if [ ! -f \${contig}_mp.txt ]; then
    NT=\$(grep -A1 "^>\${contig}\$" ${GENOME} | tail -n 1 | head -c 1)
    echo -e "\${contig}\t1\t\${NT}\t1\t0\t1\t1\t\${NT}\t0.0\t0.0\t0.0" \\
        > \${contig}_snape.txt
    exit 0
fi
ERR_FILE=\${contig}-${UNMASKED_PREFIX}.err
# Below, I redirect stderr to avoid printing many warnings about indels
snape-pooled -nchr $(($N_ADULTS*2)) -theta ${THETA} -D ${D} \\
    -priortype ${PRIOR_TYPE} -fold ${FOLD} < \${contig}_mp.txt \\
    1> \${contig}_snape.txt \\
    2> \${ERR_FILE}
status=\$?
rm \${contig}_mp.txt
# Check *.err file for anything unusual; delete if it doesn't.
NON_WARNS=\$(tail -n +5 \${ERR_FILE} | grep -v "^warning:" | wc -l)
if [ \${NON_WARNS} -ge 1 ]; then
    status=999
else
    rm \${ERR_FILE}
fi
exit \$status

EOF
chmod +x one_snape.sh




#'
#' Not sure why, but in interactive environment, replacing the below line with
#' `python3 << EOF` and removing the shebang (which directly runs the script
#' and avoids having an `all_snapes.py` file) causes errors,
#' while this version works fine. ::shrug::
#'
cat << EOF > all_snapes.py
#!/usr/bin/env python3
import subprocess as sp
import sys

def work(contig):
    """Defines the work unit on an input file"""
    cmd = "./one_snape.sh " + contig
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
    print("Starting SNAPE-pooled...")
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
    print("... Ended SNAPE-pooled")
    # Now just get the return codes:
    return_codes = [x.returncode for x in results]
    #
    for i in range(n_tasks):
        if return_codes[i] != 0:
            print("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", file=sys.stderr)
            print("SNAPE call on " + tasks[i] + " may contain errors!", file=sys.stderr)
            print("See " + tasks[i] + "-${UNMASKED_PREFIX}.err", file=sys.stderr)
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", file=sys.stderr)
    sys.exit(0)
    #
EOF

chmod +x all_snapes.py

# If it's an interactive job, uncomment lines for progress bar:
if [[ $- == *i* ]]; then sed -i 's/##> //g' all_snapes.py; fi

./all_snapes.py

rm all_snapes.py one_snape.sh





for contig in ${CONTIG_NAMES[@]}; do
    cat ${contig}_snape.txt >> ${UNMASKED_PREFIX}.txt
    rm ${contig}_snape.txt
done








#' ========================================================================
#' Convert to gSYNC file
#' Creates file ${UNMASKED_PREFIX}.sync
#' (Takes ~5 min)
python3 ./snape-files/SNAPE2SYNC.py \
    --input ${UNMASKED_PREFIX}.txt \
    --ref ${GENOME_PICKLE} \
    --output ${UNMASKED_PREFIX} \
    > ${UNMASKED_PREFIX}_2SYNC.out
check_exit_status "SNAPE2SYNC" $?


# I'm capturing output above to avoid unnecessary warnings printed to stdout
# Remove unnecessary info, print file to stdout, then delete file
grep -v "^Reference is N" ${UNMASKED_PREFIX}_2SYNC.out \
    | sed '/^[[:space:]]*$/d'
rm ${UNMASKED_PREFIX}_2SYNC.out


#' ========================================================================
#' Create masked gSYNC file
#' Creates two files:
#' - ${MASKED_FULL_PREFIX}.bed
#' - ${MASKED_FULL_PREFIX}_masked.sync
#'
#' The mask removes frequences that meet any of these:
#' - Pr(polymorphic) < 0.9
#' - read depth – too low/high
#' - within repetitive elements
#' - within 5 bp of indel
#'
#' Takes ~6 min
python3 ./snape-files/MaskSYNC_snape.py \
    --sync ${UNMASKED_PREFIX}.sync \
    --output ${MASKED_FULL_PREFIX} \
    --indel ${MP_INFO_PREFIX}.indel \
    --coverage ${MP_INFO_PREFIX}.cov \
    --te ${REPEATS} \
    --mincov ${MIN_COV} \
    --maxcov ${MAX_COV} \
    --maxsnape ${MAX_SNAPE} \
    --SNAPE
check_exit_status "MaskSYNC_snape" $?


# This removes redundant _masked ending to this file:
mv ${MASKED_FULL_PREFIX}_masked.sync ${MASKED_FULL_PREFIX}.sync

#' Create a snp file for use in npstat:
grep -v ".:.:.:.:.:." ${MASKED_FULL_PREFIX}.sync \
    | cut -f 1,2 \
    | gzip \
    > ${MASKED_FULL_PREFIX}.snp.gz



#' ========================================================================
#' Create partially masked gSYNC file for npstat
#' Differs from above bc it doesn't include filtering based on P(polymorphic)
#' Creates two files:
#' - ${MASKED_PART_PREFIX}.bed.gz
#' - ${MASKED_PART_PREFIX}_masked.sync.gz
#'
#' Takes ~5 min
python3 ./snape-files/MaskSYNC_snape.py \
    --sync ${UNMASKED_PREFIX}.sync \
    --output ${MASKED_PART_PREFIX} \
    --indel ${MP_INFO_PREFIX}.indel \
    --coverage ${MP_INFO_PREFIX}.cov \
    --te ${REPEATS} \
    --mincov ${MIN_COV} \
    --maxcov ${MAX_COV} \
    --maxsnape 0 \
    --SNAPE
check_exit_status "MaskSYNC_snape (partial)" $?


# This removes redundant _masked ending to this file:
mv ${MASKED_PART_PREFIX}_masked.sync ${MASKED_PART_PREFIX}.sync






#' ========================================================================
#' Handle output files
#'
#'

rm -r snape-files ${MP_INFO_PREFIX}* ${REPEATS} ${GENOME_PICKLE} ${GENOME}


#' Compress all output
#' Takes ~7 min
gzip *.txt *.bed *.sync


cd ..
tar -cf ${OUT_DIR}.tar ${OUT_DIR}
mv ${OUT_DIR}.tar ${PARENT_DIR_OUT}/

rm -r ${OUT_DIR}


if [ -f tany_time.sif ]; then rm tany_time.sif; fi
