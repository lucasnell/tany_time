#!/usr/bin/env python3
"""
Per-pool QC metrics from partially masked gSYNC files (one per pool, all
covering the same positions in the same order; optionally gzipped).

Output columns (tab separated, one row per pool):

  pool          pool name (input file name up to the first ".sync")
  n_sites       number of positions in the file
  n_missing     positions masked (coverage, repeat, or indel masks)
  prop_missing  n_missing / n_sites
  n_snps        positions where SNAPE P(polymorphic) >= --min-pp
  n_private     positions where this is the only input pool with
                P(polymorphic) >= --min-pp

Private SNPs depend on which pools are input, so compare them only among
pools from the same run.

usage:
pool-qc-sync.py [--min-pp 0.9] -o <output.tsv> <input sync files>

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


def snape_pp(l_split):
    """SNAPE P(polymorphic) from one split line, or -1 if masked/absent"""
    if len(l_split) < 5:
        return -1.0
    fields = l_split[4].rstrip().split(":")
    if len(fields) < 5 or fields[4] == ".":
        return -1.0
    return float(fields[4])


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Per-pool QC from syncs")
    parser.add_argument("-o", "--out", required = True,
                        help = "Output file name (.tsv)")
    parser.add_argument("--min-pp", type = float, default = 0.9,
                        help = "Minimum SNAPE P(polymorphic) for a SNP " +
                        "(default 0.9)")
    parser.add_argument("inputs", nargs = "+", help="input sync files")
    args = parser.parse_args()

    if os.path.exists(args.out):
        make_error(args.out + " already exists. Overwriting not allowed.")
    for x in args.inputs:
        if not os.path.exists(x):
            make_error(x + " not found")

    names = [os.path.basename(x).split(".sync")[0] for x in args.inputs]
    n = len(names)
    n_missing = [0] * n
    n_snps = [0] * n
    n_private = [0] * n

    print("Started (" + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ")",
          flush = True)

    files = [gzip.open(x, "rt") if x.endswith(".gz") else open(x, "r")
             for x in args.inputs]

    j = 0
    for line in it.zip_longest(*files, fillvalue = ""):
        pos0 = None
        poly = []
        for i in range(n):
            l_split = line[i].rstrip().split("\t")
            if len(l_split) < 4:
                make_error("only " + str(len(l_split)) + " columns in file " +
                           args.inputs[i] + " line " + str(j))
            if i == 0:
                pos0 = l_split[:2]
            elif l_split[:2] != pos0:
                make_error("file " + args.inputs[i] + " line " + str(j) +
                           " has a different position than " + args.inputs[0])
            if l_split[3].startswith("."):
                n_missing[i] += 1
            elif snape_pp(l_split) >= args.min_pp:
                n_snps[i] += 1
                poly.append(i)
        if len(poly) == 1:
            n_private[poly[0]] += 1
        j += 1

    for x in files:
        x.close()

    with open(args.out, "w") as out_file:
        out_file.write("pool\tn_sites\tn_missing\tprop_missing\tn_snps\t" +
                       "n_private\n")
        for i in range(n):
            out_file.write("\t".join([names[i], str(j), str(n_missing[i]),
                                      str(n_missing[i] / j if j > 0 else "NA"),
                                      str(n_snps[i]), str(n_private[i])]) +
                           "\n")

    print("Finished " + str(j) + " sites (" +
          datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ")", flush = True)
    sys.exit(0)
