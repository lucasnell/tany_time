# `04-diversity/`: crashes, immigration, and neutral diversity

ABC test of whether neutral diversity through the SN time series needs
immigration from an unsampled source, beyond local size changes that follow
the trap counts. The design (models, priors, statistics, decision rule,
sensitivity runs) is fixed in `../../_R/abc-design.md` (approved
2026-11-12). Simulations use msprime (`sim-env`, container v2.6.0);
statistics and the observed data use R (`baypass-env`, for poolfstat).

| Step | Does | Jobs | Output (`/staging/l/lnell/dna/abc/`) |
|---|---|---|---|
| `01-abc-sims` | One batch of 100 simulations for one model (`../bin/abc-sim.py` → `../bin/abc-stats.R`: pool-seq reads at the observed depths, then the 258 statistics), or the observed statistics (model `observed`) | 1,301 (`abc_info.txt`): observed; 250 reference batches and 10 held-out batches per model (M1-M5) | `sims/abc-<variant>-<model>-b<batch>.tsv.gz`; `abc-observed-stats.tsv` |

## Inputs

`abc-inputs.tar.gz`, from `_R/abc-inputs.R` (local; reads the `02-npstat`
and per-pool `01-snape` outputs of the 18 SN pools from the local data
folder), uploaded to `/staging/l/lnell/dna/abc/`. It holds pools,
trap counts, the 2,000 chosen 10-kb windows, read depth at every site of
those windows (`depth.bin`), observed second-allele counts, and the
priors (`abc-priors.csv`; `abc-priors-mu.csv` for the mutation-rate
sensitivity run).

## Checks

Local, before submitting: `../tests/abc-stats-check.R` (statistics on pools
with known diversity) and `../tests/abc-sim-check.R` (simulator and
statistics with fixed parameters; needs a python with msprime 1.4.4).
On the cluster: `./check-step.sh 04-diversity/01-abc-sims` (rows per output,
time per simulation by model, observed statistic count).

## Run

Pilot first: the observed job and batch 1 of each model (6 jobs), then the
rest:
```
cd ~/tany_time && git pull && cd _bash/04-diversity
mkdir -p logs/01-abc-sims
grep -E '^main,(observed,0|M[1-5],1),' abc_info.txt > logs/01-abc-sims/first.txt
sed 's|^queue \(.*\) from .*|queue \1 from logs/01-abc-sims/first.txt|' 01-abc-sims.sub | condor_submit -
```
