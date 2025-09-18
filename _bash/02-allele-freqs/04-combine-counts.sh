#!/bin/bash

#' This script combines all quality-passing gSYNC objects into three:
#' one for temporal samples, one for spatial samples, and one for all.
#'
#' This differs from `03-combine-snapes.sh` in that the output contains
#' *counts* of alleles by nucleotide type instead of allele frequency
#' proportion info.
#'



. /app/.bashrc
conda activate main-env


# Edit this for your system:
export SNAPE_DIR="/staging/lnell/dna/snape"


export TIME_SAMPS=(SN-1-02_S55 SN-1-07_S23 SN-1-08_S24 SN-1-09_S25 SN-1-10_S26 \
                   SN-1-13_S56 SN-1-14_S27 SN-1-97_S61 SN-1-99_S28 SN-2-00_S29 \
                   SN-2-02_S30 SN-2-07_S31 SN-2-08_S32 SN-2-15_S34 SN-2-92_S35 \
                   SN-2-97_S37 SN-2-99_S38)
export SPACE_SAMPS=(Ash-19_S5 Blik-19_S6 Blo-19_S7 Ellid-19_S8 Hrisatjorn_S11 \
                    Lei-19_S13 Lys-19_S14 Mik-19_S15 MikW-19_S16 MyBR-19_S17 \
                    MyKS-19-B_S18 MySN-19_S19 Ola-19_S20 Pernuvatn_S21 \
                    Skj-19_S22 Sva-19_S39 Tjornin_S40 Vik-19_S41)
export ALL_SAMPS=("${TIME_SAMPS[@]}" "${SPACE_SAMPS[@]}")



#' ========================================================================
#' Inputs
#' ========================================================================


for samp in ${ALL_SAMPS[@]}; do
    tar -xf ${SNAPE_DIR}/${samp}_snape.tar ${samp}_snape/${samp}_snape_masked.sync.gz \
        && mv ${samp}_snape/${samp}_snape_masked.sync.gz ./ \
        && rm -r ${samp}_snape \
        && mv ${samp}_snape_masked.sync.gz ${samp}.sync.gz
    check_exit_status "cp, rename ${samp}" $?
done




#' Combine multiple gSYNC files (optionally gzipped) into one.
#' It assumes they're all in the same order, as they should be because this
#' format should include a line per bp in the reference.
#' It also removes any instances of a 5th column in the output, which
#' can happen when converting from SNAPE-pooled output to gSYNC.
#'
#' For use in `poolfstat`, you should also run this to remove any missing data:
#' gunzip -c POOLS.sync.gz | grep -v ".:.:.:.:.:." | gzip > POOLS_noblanks.sync.gz
#'
#' usage:
#' ./combine-counts.py \
#'     -o <output file name> \
#'     <input sync files as *.sync or *.sync.gz>

cat << EOF > combine-counts.py
#!/usr/bin/env python3

import sys
import os.path
import gzip
import argparse
import itertools as it
from datetime import datetime


def make_error(err_msg):
    sys.stderr.write("ERROR: " + err_msg + "\n")
    sys.exit(1)



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Combine sync files")
    parser.add_argument("-o", "--out", required = True,
                        help = "Output file name")
    parser.add_argument("inputs", nargs = "+", help="input sync files")

    args = parser.parse_args()
    out = args.out
    inputs = args.inputs

    if not out.endswith((".sync", ".sync.gz")):
        make_error("Strange suffix to output file name. Exiting.")
    out_dir_file = os.path.split(out)
    if out_dir_file[0] != "" and not os.path.exists(out_dir_file[0]):
        make_error("output dir (" + out_dir_file[0] + ") not found")
    if os.path.exists(out):
        make_error(out + " already exists. Overwriting not allowed.\n")

    if len(inputs) < 2:
        make_error("At least 2 input files required")

    gzipped = all(x.endswith(".sync.gz") for x in inputs)
    if not (gzipped or all(x.endswith(".sync") for x in inputs)):
        make_error("All inputs must have the same extension, and it must " +
                   "be *.sync.gz or *.sync")

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print("Combining sync files into " + out)

    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("started reading/writing (" + time_str + ")")

    if gzipped:
        files = [gzip.open(x, "rt") for x in inputs]
    else:
        files = [open(x, "r") for x in inputs]

    if out.endswith(".gz"):
        out_file = gzip.open(out, "wt")
    else:
        out_file = open(out, "w")

    j = 0
    for line in it.zip_longest(*files, fillvalue=""):
        l_split = line[0].rstrip().split("\t")
        if len(l_split) < 4:
            make_error("only " + str(len(l_split)) + " columns in file " +
                       inputs[0] + " line " + str(j))
        zline = l_split[:4]
        for i in range(1, len(line)):
            l_split = line[i].split("\t")
            if len(l_split) < 4:
                make_error("only " + str(len(l_split)) + " columns in file " +
                           inputs[0] + " line " + str(j))
            zline.append(l_split[3].rstrip())
        new_line = "\t".join(zline) + "\n"
        b = out_file.write(new_line)
        j += 1

    for x in files:
        x.close()
    out_file.close()

    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("... finished (" + time_str + ")")

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    sys.exit(0)

EOF

chmod +x combine-counts.py



#' ========================================================================
#' Combine files
#' ========================================================================


for type in TIME SPACE ALL; do

    #' Using indirect reference (notice `!`) to retrieve one of the
    #' *_SAMPS arrays:
    SAMPS="${type}_SAMPS[@]"
    SAMPS=("${!SAMPS}")

    # Names of output files:
    OUT_BASE=count_masked_$(echo "${type}" | tr '[:upper:]' '[:lower:]')
    OUT=${OUT_BASE}.sync.gz
    OUT_NAMES=${OUT_BASE}.names.gz
    OUT_NOBLANKS=${OUT_BASE}_noblanks.sync.gz

    # Combine snape files (takes ~17 min)
    ./combine-counts.py -o ${OUT%.gz} ${SAMPS[@]/%/.sync.gz}
    check_exit_status "combine-counts.py ($type)" $?

    # Create file of names
    for samp in "${SAMPS[@]}"; do
        echo ${samp} >> ${OUT_NAMES%.gz}
    done

    # Create count file with no blanks
    # Note that this version removes lines with 6 '.', not 7 like in
    # `combine-snapes.sh`
    cat ${OUT%.gz} | grep -Fv ".:.:.:.:.:." | gzip > ${OUT_NOBLANKS}
    check_exit_status "make_no_blanks ($type)" $?

    # compress output
    gzip ${OUT_NAMES%.gz} \
        && gzip ${OUT%.gz}
    check_exit_status "gzip out and names ($type)" $?

    # create output tar file
    # Note that I'm copying (not moving) to staging bc this file will be
    # transferred to ResearchDrive upon exit
    mkdir ${OUT_BASE}
    mv ${OUT} ${OUT_NAMES} ${OUT_NOBLANKS} ./${OUT_BASE}/ \
        && tar -cf ${OUT_BASE}.tar ${OUT_BASE} \
        && cp ${OUT_BASE}.tar ${SNAPE_DIR}/
    check_exit_status "tar, copy over output files ($type)" $?

    # remove unnecessary folder of outputs
    rm -r ${OUT_BASE}

done




rm *.sync.gz *.py
if [ -f tany_time.sif ]; then rm tany_time.sif; fi

exit 0
