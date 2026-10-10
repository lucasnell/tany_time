#!/usr/bin/env python3
"""
Combine multiple gSYNC files (optionally gzipped) into one.

It assumes they're all in the same order, as they should be because this
format includes a line per bp in the reference.

By default, the output contains allele *counts* (4th column of each input),
and any 5th column in the inputs is dropped. With `--snape`, the output
instead contains the SNAPE-pooled info (5th column of each input);
inputs lacking a 5th column at a site get ".:.:.:.:.:.:." there.

With `--min-pp P`, a site is written only if at least one input has
SNAPE P(polymorphic) >= P there (5th field of the SNAPE info column);
inputs without SNAPE info at a site (i.e., masked) don't count. Output then
no longer has one line per bp.

To remove sites with missing data in any pool:
gunzip -c POOLS.sync.gz | grep -Fv ".:.:.:.:.:." | gzip > POOLS_noblanks.sync.gz
(This pattern also matches the 7-field SNAPE placeholder.)

usage:
combine-syncs.py [--snape] [--min-pp P] -o <output file name> <input sync files>

"""

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
    parser.add_argument("--snape", action = "store_true",
                        help = "Combine SNAPE info (5th column) instead " +
                        "of allele counts (4th column)")
    parser.add_argument("--min-pp", type = float, default = None,
                        help = "Keep only sites where at least one input " +
                        "has SNAPE P(polymorphic) >= this value")
    parser.add_argument("inputs", nargs = "+", help="input sync files")

    args = parser.parse_args()
    out = args.out
    inputs = args.inputs
    snape = args.snape
    min_pp = args.min_pp

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
    print("Combining sync files (" + ("SNAPE info" if snape else "counts") +
          ") into " + out)

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

    def get_field(l_split):
        """Column to keep from one input line (already split by tab)"""
        if not snape:
            return l_split[3].rstrip()
        if len(l_split) < 5:
            return ".:.:.:.:.:.:."
        return l_split[4].rstrip()

    def max_pp(l_split):
        """SNAPE P(polymorphic) for one input line, or -1 if absent"""
        if len(l_split) < 5:
            return -1.0
        fields = l_split[4].rstrip().split(":")
        if len(fields) < 5 or fields[4] == ".":
            return -1.0
        return float(fields[4])

    j = 0
    n_written = 0
    for line in it.zip_longest(*files, fillvalue=""):
        zline = []
        pp = -1.0
        for i in range(len(line)):
            l_split = line[i].rstrip().split("\t")
            if len(l_split) < 4:
                make_error("only " + str(len(l_split)) + " columns in file " +
                           inputs[i] + " line " + str(j))
            if i == 0:
                zline = l_split[:3]
            elif l_split[:2] != zline[:2]:
                make_error("file " + inputs[i] + " line " + str(j) +
                           " has a different position than " + inputs[0])
            zline.append(get_field(l_split))
            if min_pp is not None:
                pp = max(pp, max_pp(l_split))
        j += 1
        if min_pp is not None and pp < min_pp:
            continue
        new_line = "\t".join(zline) + "\n"
        b = out_file.write(new_line)
        n_written += 1

    for x in files:
        x.close()
    out_file.close()

    print("wrote " + str(n_written) + " of " + str(j) + " sites")
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("... finished (" + time_str + ")")

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    sys.exit(0)
