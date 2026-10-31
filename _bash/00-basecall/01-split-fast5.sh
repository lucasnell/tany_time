#!/bin/bash

#' Split the raw nanopore signal archive (`ives_fast5.tar.bz2`, 330 GB;
#' ~800 multi-read FAST5 files of ~400 MB) into uncompressed tar chunks of
#' N_PER_CHUNK FAST5 files, one per `02-dorado` GPU job. Checks the
#' archive's MD5 in the same pass (expected value from the README next to
#' the archive on the external drive) and fails (exit 3) on a mismatch.
#' Single-threaded bz2 decompression; expect several hours.
#'
#' Output (`/staging/l/lnell/dna/ont/fast5_chunks/`): `chunk_NNN.tar` and
#' `manifest.tsv` (chunk, FAST5 file).
#'


. /app/.bashrc
. bin/helpers.sh || exit 1
conda activate main-env

# Edit these for your system:
export ARCHIVE="/staging/l/lnell/dna/ont/ives_fast5.tar.bz2"
export PARENT_DIR_OUT="/staging/l/lnell/dna/ont/fast5_chunks"

export N_PER_CHUNK=20
export EXPECTED_MD5="8dcd7ef5bac44a8f9146bceef059f29f"

require_files ${ARCHIVE}
mkdir -p ${PARENT_DIR_OUT}

python3 -u << 'PYEOF'
import bz2, hashlib, os, shutil, sys, tarfile

arch = os.environ["ARCHIVE"]
out = os.environ["PARENT_DIR_OUT"]
n_per = int(os.environ["N_PER_CHUNK"])
expected = os.environ["EXPECTED_MD5"]


class Hashing:
    """File wrapper that updates an MD5 with every byte read."""
    def __init__(self, f):
        self.f = f
        self.h = hashlib.md5()
    def read(self, n=-1):
        b = self.f.read(n)
        self.h.update(b)
        return b


raw = Hashing(open(arch, "rb"))
stream = tarfile.open(fileobj=bz2.BZ2File(raw), mode="r|")
os.makedirs("f5", exist_ok=True)
manifest = open("manifest.tsv", "w")
chunk, members, n = 0, [], 0


def flush():
    global chunk, members
    if not members:
        return
    name = "chunk_%03d.tar" % chunk
    with tarfile.open(name, "w") as ct:
        for m in members:
            ct.add(m, arcname=os.path.basename(m))
    shutil.move(name, os.path.join(out, name))
    for m in members:
        manifest.write("%s\t%s\n" % (name, os.path.basename(m)))
        os.remove(m)
    manifest.flush()
    print("wrote", name, len(members), "files", flush=True)
    chunk += 1
    members = []


for m in stream:
    base = os.path.basename(m.name)
    # skip non-FAST5 entries and macOS metadata files ("._*")
    if not m.isfile() or not base.endswith(".fast5") or base.startswith("._"):
        continue
    dst = os.path.join("f5", base)
    with open(dst, "wb") as o:
        shutil.copyfileobj(stream.extractfile(m), o)
    members.append(dst)
    n += 1
    if len(members) == n_per:
        flush()
flush()
manifest.close()

# read any bytes left after the tar stream so the MD5 covers the whole file
while raw.read(1 << 24):
    pass
md5 = raw.h.hexdigest()
print("fast5 files:", n, "chunks:", chunk)
print("archive MD5:", md5, "expected:", expected)
sys.exit(0 if md5 == expected else 3)
PYEOF
check_exit_status "split archive (exit 3 = MD5 mismatch)" $?

mv manifest.tsv ${PARENT_DIR_OUT}/
check_exit_status "move manifest to staging" $?

safe_exit 0
