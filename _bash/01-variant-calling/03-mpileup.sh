#!/bin/bash

#' This script does the following:
#' - mark and remove duplicates using `Picard MarkDuplicates`
#' - re-align sequences in the proximity of indels with `IndelRealigner` and
#'   `RealignerTargetCreator` in `GATK`
#' - `samtools mpileup`



. /app/.bashrc
conda activate main-env

export THREADS=$(count_threads)

# Edit these for your system:
export PARENT_DIR_IN="/staging/lnell/dna/bwa"
export PARENT_DIR_OUT="/staging/lnell/dna/mpileup"
export GENOME_FULL_PATH="/staging/lnell/Tgraci_assembly.fasta.gz"





#' ========================================================================
#' Inputs
#' ========================================================================

export READ_BASE=$1

export IN_BAM=${READ_BASE}_bwa.bam
export GENOME=$(basename ${GENOME_FULL_PATH%.gz})

if [ ! -f ${PARENT_DIR_IN}/${IN_BAM} ]; then
    echo "${PARENT_DIR_IN}/${IN_BAM} does not exist! " 1>&2
    safe_exit 111
fi
if [ ! -f ${GENOME_FULL_PATH} ]; then
    echo "${GENOME_FULL_PATH} does not exist! " 1>&2
    safe_exit 222
fi


#' ========================================================================
#' Window pileup script
#' ========================================================================

cat << EOF > window-pileup.py
#!/usr/bin/env python3
"""
Summarize coverage from mpileup using a sliding window.

The mpileup file can optionally only contain 3 columns: sequence name,
position, and coverage. This can be a useful format for reducing file sizes.
This script will detect which format applies based on the first line.

usage:
./window-mpileup.py \
    -r <reference assembly> \
    -s <step size (default 50)> \
    -w <window size (default 100)> \
    <input mpileup as *.txt or *.txt.gz>

"""



import sys
import os.path
import gzip
import argparse
import numpy as np
import math
import statistics as stats
from datetime import datetime



def make_error(err_msg):
    sys.stderr.write("ERROR: " + err_msg + "\n")
    sys.exit(1)



# read first line of mpileup file to see if it's an okay format
# and whether it has 3 or >= 5 columns
# returns whether it has 3 columns
def test_mpileup(in_mpileup):

    if in_mpileup.endswith(".gz"):
        with gzip.open(in_mpileup,"rt") as mp_file:
            first_line = mp_file.readline().rstrip().split("\t")
    else:
        with open(in_mpileup, "r") as mp_file:
            first_line = mp_file.readline().rstrip().split("\t")

    if len(first_line) != 3 and len(first_line) < 5:
        make_error("mpileup file must have 3 or >= 5 columns")

    try:
        tmp = int(first_line[1])
    except ValueError:
        make_error("the mpileup file's 2nd col must be an integer")
    except:
        make_error("Something else went wrong")

    three_cols = len(first_line) == 3

    if three_cols:
        try:
            tmp = int(first_line[2])
        except ValueError:
            make_error("For 3-column mpileup file, 3rd col must be an integer")
        except:
            make_error("Something else went wrong")
    else:
        try:
            tmp = int(first_line[3])
        except ValueError:
            make_error("For mpileup file with > 3 cols, 4th col must be an integer")
        except:
            make_error("Something else went wrong")

    return three_cols



def summarize_ref(reference):

    if reference.endswith(".gz"):
        ref_file = gzip.open(reference,"rt")
    else:
        ref_file = open(reference, "r")

    sizes = []
    names = []
    i = -1

    for line in ref_file:
        if line.startswith(">"):
            sizes.append(0)
            seq_name = line[1:]
            seq_name = seq_name.strip(" \t\n\r")
            names.append(seq_name)
            i += 1
        else:
            if len(sizes) == 0:
                make_error("FASTA doesn't start with header. Exiting.\n")
            # we don't want to include newline at the end in our counts:
            llen = len(line) - 1
            sizes[i] += llen

    ref_file.close()

    return sizes, names



def one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends, names):
    """
    Summarize one window. Updates the output file directly.
    """

    if len(cov_ij) == 0:
        out_line = names[seq_i] + "\t"
        out_line += str(starts[seq_i][win_j]) + "\t"
        out_line += str(ends[seq_i][win_j]) + "\t"
        out_line += "0\t" + "0\t" + "0.0\t" + "0.0\n"
        b = out_file.write(out_line)
        return

    size_ij = np.int64(ends[seq_i][win_j] - starts[seq_i][win_j] + 1)
    if len(cov_ij) > size_ij:
        make_error("For window " + str(win_j+1) +
                   " in sequence '" + names[seq_i] +
                   "', the number of coverage values exceeds " +
                   "the window size!\n")

    n_zeros = size_ij - np.int64(len(cov_ij))
    for i in range(n_zeros):
        cov_ij.append(np.float64(0))

    min_ij = np.int64(min(cov_ij))
    max_ij = np.int64(max(cov_ij))
    mean_ij = stats.mean(cov_ij)
    median_ij = stats.median(cov_ij)

    out_line = names[seq_i] + "\t"
    out_line += str(starts[seq_i][win_j]) + "\t"
    out_line += str(ends[seq_i][win_j]) + "\t"
    out_line += str(min_ij) + "\t"
    out_line += str(max_ij) + "\t"
    out_line += str(mean_ij) + "\t"
    out_line += str(median_ij) + "\n"
    b = out_file.write(out_line)

    return


# For testing:
# in_mpileup = "Lys-19_S14_bwa_mpileup.txt.gz"
# reference = "tany_scaffolds.fasta.gz"
# step = 50
# window = 100



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Summarize mpileup")
    parser.add_argument("-r", "--reference", required = True,
                        help = "reference assembly file name")
    parser.add_argument("-s", "--step", type=int, default=50,
                        help = "sliding window step size (default = 50)")
    parser.add_argument("-w", "--window", type=int, default=100,
                        help = "max coverage before binning (default = 100)")
    parser.add_argument("-o", "--out_fn", required = False,
                        help = "output file name (defaults to "+
                        "${INPUT_NAME/.txt/_window.txt})")
    parser.add_argument("in_mpileup", help="input mpileup file name")

    args = parser.parse_args()
    in_mpileup = args.in_mpileup
    reference = args.reference
    step = args.step
    window = args.window


    # ---------------
    # Check files:
    # ---------------
    if not os.path.exists(in_mpileup):
        make_error(in_mpileup + " not found\n")
    mp_suffs = (".txt.gz", ".txt")
    if not in_mpileup.endswith(mp_suffs):
        make_error("Only *.txt or *.txt.gz extension allowed. Exiting.\n")
    fasta_suffs = (".fasta.gz", ".fasta", ".fa.gz", ".fa")
    if not reference.endswith(fasta_suffs):
        make_error("Strange suffix to input FASTA file. Exiting.\n")
    if args.out_fn:
        if not args.out_fn.endswith((".txt", ".txt.gz")):
            make_error("Strange suffix to output file name. Exiting.")
        out_dir_file = os.path.split(args.out_fn)
        if out_dir_file[0] != "" and not os.path.exists(out_dir_file[0]):
            make_error("output dir (" + out_dir_file[0] + ") not found")
        out_fn = args.out_fn
    else:
        out_fn = in_mpileup.replace(".txt", "_window.txt")
    if os.path.exists(out_fn):
        make_error(out_fn + " already exists. Overwriting not allowed.\n")

    three_cols = test_mpileup(in_mpileup)


    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print("Sliding windows of coverage for " + in_mpileup + "\n")

    # --------------
    # Read reference file:
    # --------------
    sizes, names = summarize_ref(reference)
    num_seqs = len(sizes)

    # --------------
    # Make output objects:
    # --------------
    # vector of windows per sequence:
    win_vec = np.array([math.ceil((ss - window + 1) / step) for ss in sizes])

    # starting and ending points (inclusive) for each window:
    starts = [np.arange(0, ss-window+1, step) for ss in sizes]
    ends = [np.append(np.arange(window-1, ss-step, step), ss-1) for ss in sizes]



    if in_mpileup.endswith(".gz"):
        mp_file = gzip.open(in_mpileup,"rt")
    else:
        mp_file = open(in_mpileup, "r")
    if out_fn.endswith(".gz"):
        out_file = gzip.open(out_fn,"wt")
    else:
        out_file = open(out_fn, "w")

    out_line = "scaff\t" + "start\t" + "end\t" + "min\t" + "max\t"
    out_line += "mean\t" + "median\n"
    # Assigning to `b` so it doesn't print # bytes printed
    b = out_file.write(out_line)

    # indices for sequence and window
    seq_i = np.int64(0)
    win_j = np.int64(0)
    # to calc. median, we need to retain all coverage values within a window:
    cov_ij = []

    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("starting processing files (" + time_str + ")")

    for line in mp_file:

        if three_cols:
            scaff, pos, cov = line.split("\t")
        else:
            split_line = line.split("\t")
            scaff = split_line[0]
            pos = split_line[1]
            cov = split_line[3]
        pos = np.int64(pos)
        # mpileup uses 1-based indices
        pos -= 1
        cov = np.float64(cov)

        if scaff != names[seq_i]:
            """
            If we need to move to a new sequence, then summarize this window
            and, if applicable, any remaining in the sequence.
            Then do this for all windows in the sequences between this one
            and the one you need.
            """
            one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends, names)
            cov_ij = []
            win_j += 1
            while scaff != names[seq_i]:
                while win_j < len(ends[seq_i]):
                    one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends,
                                  names)
                    win_j += 1
                seq_i += 1
                win_j = np.int64(0)
                if seq_i >= len(names):
                    make_error("assembly and mpileup sequence names " +
                               "do not match or are not in the same order!")

        if pos > ends[seq_i][win_j]:
            """
            If we need to change window within the sequence, then similarly
            summarize windows until you get to the one you want.
            """
            one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends, names)
            cov_ij = []
            win_j += 1
            while pos > ends[seq_i][win_j]:
                one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends,
                                  names)
                win_j += 1
                if win_j >= len(ends[seq_i]):
                    make_error("position " + str(pos) +
                               " exceeds highest possible index for " +
                               "sequence '" + names[seq_i] + "' (" +
                               str(sizes[seq_i]-1)  + ")\n")

        # Once you're at the sequence and window you need, then append this
        # coverage to `cov_ij`
        cov_ij.append(cov)

    # Summarize last window
    one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends, names)
    # Summarize any remaining empty windows:
    cov_ij = []
    win_j += 1
    while seq_i < len(names):
        while win_j < len(ends[seq_i]):
            one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends, names)
            win_j += 1
        seq_i += 1
        win_j = np.int64(0)

    mp_file.close()
    out_file.close()

    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("... finished processing files (" + time_str + ")")

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    sys.exit(0)

EOF

chmod +x window-pileup.py




#' ========================================================================
#' Outputs
#' ========================================================================


# Final files / directories
export OUT_DIR=${READ_BASE}_mpileup
export OUT_FILE=${READ_BASE}_mpileup.txt.gz
# Intermediates:
export MARKDUP_OUT=${IN_BAM/.bam/_nodups.bam}
export REALIGNED_OUT=${MARKDUP_OUT/.bam/_realigned.bam}

# if [ -f ${PARENT_DIR_OUT}/${OUT_DIR}.tar.gz ] && [ -f ${PARENT_DIR_OUT}/${OUT_FILE} ]; then
#     echo "Output files already exist" 1>&2
#     safe_exit 0
# fi


#' ========================================================================
#' Prep for downstream steps.
#' ========================================================================

mkdir ${OUT_DIR}
cd ${OUT_DIR}
mv ../window-pileup.py ./

cp ${PARENT_DIR_IN}/${IN_BAM} ./
check_exit_status "move BAM" $?

cp ${GENOME_FULL_PATH} ./ && \
    gunzip ${GENOME}.gz
check_exit_status "move, gunzip genome" $?


# Needed downstream
samtools faidx --length 80 ${GENOME}



#' ========================================================================
#' Mark and remove duplicates
#' ========================================================================

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
conda activate gatk-env

gatk CreateSequenceDictionary R=${GENOME}
check_exit_status "Picard_CreateSequenceDictionary" $?

gatk MarkDuplicates \
    REMOVE_DUPLICATES=true \
    I=${IN_BAM} \
    O=${MARKDUP_OUT} \
    M=${MARKDUP_OUT/.bam/_report.txt} \
    VALIDATION_STRINGENCY=SILENT \
    VERBOSITY=WARNING
check_exit_status "Picard_MarkDuplicates" $?

conda deactivate
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

rm ${IN_BAM}

call_bam_stats ${MARKDUP_OUT} "(nodups)"

samtools index -@ $(($THREADS - 1)) ${MARKDUP_OUT}
check_exit_status "samtools index (nodups)" $?


#' ========================================================================
#' Realign around indels
#' ========================================================================

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
conda activate gatk3-env

GenomeAnalysisTK \
    -T RealignerTargetCreator \
    -nt ${THREADS} \
    -R ${GENOME} \
    -I ${MARKDUP_OUT} \
    -o ${MARKDUP_OUT/.bam/.intervals} \
    -log ${MARKDUP_OUT/.bam/.intervals.log}
check_exit_status "RealignerTargetCreator" $?

# This part can take a while (9+ hours)
GenomeAnalysisTK \
    -T IndelRealigner \
    -R ${GENOME} \
    -I ${MARKDUP_OUT} \
    -targetIntervals ${MARKDUP_OUT/.bam/.intervals} \
    -o ${REALIGNED_OUT}
check_exit_status "IndelRealigner" $?

conda deactivate
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

rm ${MARKDUP_OUT}*

call_bam_stats ${REALIGNED_OUT} "(realigned)"


#' ========================================================================
#' mpileup
#' ========================================================================
samtools mpileup -B -f ${GENOME} ${REALIGNED_OUT} > ${OUT_FILE%.gz}
check_exit_status "samtools mpileup" $?

gzip ${OUT_FILE%.gz}

# Sliding window (500bp with step size of 250bp) of coverage for plotting
./window-mpileup.py -r ${GENOME} -s 250 -w 500 ${OUT_FILE}
check_exit_status "window-mpileup.py" $?
rm window-pileup.py

# Change default naming to specify window size:
mv ${OUT_FILE/.txt/_window.txt} ${OUT_FILE/.txt/_win500.txt}

rm ${GENOME/.fasta/}*




#' ========================================================================
#' Handle output files
#' ========================================================================

#'
#' Main mpileup file is separated from directory containing logs.
#' Both the mpileup file and the archived directory will be copied to
#' both staging and (upon exit) ResearchDrive.
#'
mv ${OUT_FILE} ../
cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

# Copying, not moving so that they will be moved to ResearchDrive after exit
cp ${OUT_DIR}.tar.gz ${PARENT_DIR_OUT}/
cp ${OUT_FILE} ${PARENT_DIR_OUT}/

rm -r ${OUT_DIR}

#' Do NOT use `safe_exit` because that'll remove the files before they
#' can be moved to ResearchDrive.
exit 0
