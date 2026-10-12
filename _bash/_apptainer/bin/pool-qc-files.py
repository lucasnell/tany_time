#!/usr/bin/env python3
"""
Per-pool QC metrics from upstream pipeline outputs. Prints one tab-separated
row (with --header, prints the header instead). Missing files give "NA".

Columns:
  pool             pool name
  fastp_gc         GC fraction of reads after fastp filtering (fastp.json)
  fastp_dup        fastp's duplication rate estimate (fastp.json)
  picard_dup       PERCENT_DUPLICATION from Picard MarkDuplicates (fraction)
  map_rate         mapped / total records in the unfiltered alignments
                   (bamtools stats of all_*_bwa_{merged,unmerged}; merged
                   pairs count as one record; secondary alignments from
                   `bwa mem -M` count in both)
  usable_rate      records in the final filtered BAM (MAPQ >= 20, primary)
                   / total records in the unfiltered alignments
  target_rate      records in the final filtered BAM on T. gracilentus
                   contigs / total records in the unfiltered alignments
                   (from `<pool>_species-counts.tsv`; competitive mapping)
  nontarget_rate   same, for all other species in the combined reference

usage:
pool-qc-files.py --header
pool-qc-files.py --pool P --fastp fastp.json --bwa-dir <dir with .stats> \\
    --markdup <MarkDuplicates report>

"""

import sys
import os.path
import json
import argparse


def fmt(x):
    return "NA" if x is None else str(x)


def fastp_stats(fn):
    if fn is None or not os.path.exists(fn):
        return None, None
    with open(fn) as f:
        d = json.load(f)
    gc = d.get("summary", {}).get("after_filtering", {}).get("gc_content")
    dup = d.get("duplication", {}).get("rate")
    return gc, dup


def bamtools_counts(fn):
    """(total, mapped) records from a `bamtools stats` output file"""
    if not os.path.exists(fn):
        return None, None
    total = mapped = None
    with open(fn) as f:
        for line in f:
            if line.startswith("Total reads:"):
                total = int(line.split(":")[1].split()[0])
            elif line.startswith("Mapped reads:"):
                mapped = int(line.split(":")[1].split()[0])
    return total, mapped


def species_counts(fn):
    """dict of species -> reads from `<pool>_species-counts.tsv`"""
    if not os.path.exists(fn):
        return None
    out = {}
    with open(fn) as f:
        next(f)
        for line in f:
            sp, n = line.rstrip().split("\t")
            out[sp] = int(n)
    return out


def picard_dup(fn):
    if fn is None or not os.path.exists(fn):
        return None
    with open(fn) as f:
        lines = f.read().splitlines()
    for k, line in enumerate(lines):
        if line.startswith("## METRICS CLASS"):
            header = lines[k + 1].split("\t")
            values = lines[k + 2].split("\t")
            if "PERCENT_DUPLICATION" in header:
                return float(values[header.index("PERCENT_DUPLICATION")])
    return None


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Per-pool QC from files")
    parser.add_argument("--header", action = "store_true")
    parser.add_argument("--pool")
    parser.add_argument("--fastp")
    parser.add_argument("--bwa-dir")
    parser.add_argument("--markdup")
    args = parser.parse_args()

    cols = ["pool", "fastp_gc", "fastp_dup", "picard_dup", "map_rate",
            "usable_rate", "target_rate", "nontarget_rate"]
    if args.header:
        print("\t".join(cols))
        sys.exit(0)
    if args.pool is None:
        sys.stderr.write("ERROR: --pool is required\n")
        sys.exit(1)

    gc, fdup = fastp_stats(args.fastp)
    pdup = picard_dup(args.markdup)

    map_rate = usable_rate = target_rate = nontarget_rate = None
    if args.bwa_dir is not None:
        p = os.path.join(args.bwa_dir, args.pool)
        tm = [bamtools_counts(os.path.join(args.bwa_dir, "all_" + args.pool +
                                           "_bwa_" + x + ".stats"))
              for x in ("merged", "unmerged")]
        if all(t is not None and m is not None for t, m in tm):
            total = sum(t for t, m in tm)
            mapped = sum(m for t, m in tm)
            if total > 0:
                map_rate = mapped / total
                final_total, _ = bamtools_counts(p + "_bwa.stats")
                if final_total is not None:
                    usable_rate = final_total / total
                sp = species_counts(p + "_species-counts.tsv")
                if sp is not None:
                    target_rate = sp.get("target", 0) / total
                    nontarget_rate = sum(v for k, v in sp.items()
                                         if k != "target") / total

    print("\t".join([args.pool] + [fmt(x) for x in (gc, fdup, pdup, map_rate,
                                                     usable_rate, target_rate,
                                                     nontarget_rate)]))
    sys.exit(0)
