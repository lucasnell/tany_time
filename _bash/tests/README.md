# `tests/`: local checks of pipeline behavior

Scripts run locally (not on the cluster) from the repo root. Each writes
its result to `_data/`.

| Script | Checks | Output |
|---|---|---|
| `npstat-maxcov.sh` | npstat output is identical when `-maxcov` is lowered to the maximum coverage present, and how much faster that is (`02-allele-freqs/02-npstat` relies on this) | `_data/npstat-maxcov-check.txt` |
