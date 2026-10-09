#!/usr/bin/env python3
"""
Run a shell command once per item (e.g., per contig) in parallel.

`{}` in the command is replaced by each item. Return codes that aren't 0
are listed on stderr. Exits with 0 if all commands succeeded and 1
otherwise, so callers decide whether a failure is fatal.

usage:
run-parallel.py -t <threads> -c '<command with {}>' <item1> <item2> ...

example:
run-parallel.py -t 8 -c './one_snape.sh {}' CONTIG1 CONTIG2 CONTIG3

"""

import sys
import argparse
import subprocess as sp
from concurrent.futures import ThreadPoolExecutor, as_completed


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Run command per item")
    parser.add_argument("-t", "--threads", type = int, required = True,
                        help = "number of commands to run at once")
    parser.add_argument("-c", "--command", required = True,
                        help = "shell command, with {} replaced by each item")
    parser.add_argument("items", nargs = "+", help = "items to iterate over")
    args = parser.parse_args()

    if args.threads < 1:
        sys.stderr.write("ERROR: --threads must be >= 1\n")
        sys.exit(1)
    if "{}" not in args.command:
        sys.stderr.write("ERROR: --command must contain {}\n")
        sys.exit(1)

    n = len(args.items)
    print("Starting " + str(n) + " jobs with " + str(args.threads) +
          " at a time...", flush = True)

    def work(item):
        return sp.run(args.command.replace("{}", item), shell = True).returncode

    codes = {}
    # Threads (not processes) are enough because each one only waits on
    # its own subprocess.
    with ThreadPoolExecutor(max_workers = args.threads) as ex:
        futures = {ex.submit(work, x): x for x in args.items}
        for k, fut in enumerate(as_completed(futures)):
            item = futures[fut]
            codes[item] = fut.result()
            print("[" + str(k + 1) + "/" + str(n) + "] " + item +
                  " (return code " + str(codes[item]) + ")", flush = True)

    failed = [x for x in args.items if codes[x] != 0]
    for x in failed:
        sys.stderr.write("Command failed for " + x + " (return code " +
                         str(codes[x]) + ")\n")
    print("... finished", flush = True)

    sys.exit(1 if len(failed) > 0 else 0)
