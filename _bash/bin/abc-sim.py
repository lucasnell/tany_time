#!/usr/bin/env python3
"""
One ABC simulation (design: `_R/abc-design.md`): draws parameters for a
model from the priors, simulates the SN time series with msprime, and
writes derived-allele counts per pool at segregating sites.

Time runs backwards from 2015b (msprime time 0). Trap generation g
(0 = 1977a, 77 = 2015b) is msprime time 77 - g; SN's size between times
k and k + 1 is Ne of generation 77 - k, and Ne_anc from time 78 back.
Discrete-time Wright-Fisher in generations where SN's size is < 1e4, the
standard coalescent otherwise (see "ancestry model" below).

Loci: `--n-loci` independent loci of `--locus-len` bp, simulated
separately, matching the observed 10-kb windows one to one.

Outputs (`--out` prefix):
  <out>.params.tsv   parameter values (name, value)
  <out>.pos.bin      int32, position = locus * locus_len + offset
  <out>.counts.bin   int16, sites x pools (row-major), derived counts
"""

import argparse
import math

import msprime
import numpy as np


N_GEN = 78  # 1977a to 2015b
DTWF_MAX = 1e4  # sizes below this use discrete-time Wright-Fisher
GHOST_DTWF_MAX = 1e5  # source size cap in Wright-Fisher generations

# Parameters each model draws (scale and floor for crash-driven local size;
# Ne_c for constant; Ne_s, T_s and m or I for immigration):
MODEL_PARS = {
    "M1": ["scale", "floor", "Ne_anc"],
    "M2": ["scale", "floor", "Ne_anc", "Ne_s", "T_s", "m"],
    "M3": ["Ne_c", "Ne_anc", "Ne_s", "T_s", "m"],
    "M4": ["Ne_c", "Ne_anc"],
    "M5": ["scale", "floor", "Ne_anc", "Ne_s", "T_s", "I"],
}


def read_priors(fn):
    """CSV with columns param,lo,hi (log-uniform bounds)."""
    pr = {}
    with open(fn) as f:
        header = f.readline().strip().split(",")
        for line in f:
            x = dict(zip(header, line.strip().split(",")))
            pr[x["param"]] = (float(x["lo"]), float(x["hi"]))
    return pr


def draw(rng, lo, hi, widen):
    lo, hi = lo / widen, hi * widen
    return float(math.exp(rng.uniform(math.log(lo), math.log(hi))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, choices=sorted(MODEL_PARS))
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--priors", required=True)
    ap.add_argument("--trap", required=True,
                    help="trap counts, one per generation (78 lines, 1977a first)")
    ap.add_argument("--pools", required=True,
                    help="TSV: pool, gen (0-77), n_adults; in output order")
    ap.add_argument("--n-loci", type=int, required=True)
    ap.add_argument("--locus-len", type=int, default=10000)
    ap.add_argument("--mu", type=float, default=2.1e-9)
    ap.add_argument("--rec", type=float, default=1e-8)
    # Sensitivity runs (`_R/abc-design.md`):
    ap.add_argument("--ne-form", choices=["linear", "power"], default="linear")
    ap.add_argument("--widen", type=float, default=1.0,
                    help="divide lower and multiply upper prior bounds")
    ap.add_argument("--older", action="store_true",
                    help="one extra size change before 1977")
    # Checks and posterior-predictive runs:
    ap.add_argument("--fix", action="append", default=[],
                    help="name=value: set a parameter instead of drawing it")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    priors = read_priors(args.priors)

    trap = np.loadtxt(args.trap)
    assert trap.shape == (N_GEN,)
    pools = []
    with open(args.pools) as f:
        f.readline()
        for line in f:
            p, g, n = line.strip().split("\t")
            pools.append((p, int(g), int(n)))

    # ---- parameters ----
    pars = list(MODEL_PARS[args.model])
    crash = "scale" in pars
    if crash and args.ne_form == "power":
        pars = [x for x in pars if x != "floor"] + ["b"]
    if args.older:
        pars += ["Ne_old", "T_old"]
    th = {x: draw(rng, *priors[x], args.widen) for x in pars}
    if "b" in th:
        # b is bounded by (0.1, 1) and isn't widened:
        th["b"] = draw(rng, *priors["b"], 1.0)
    if "m" in th:
        # A proportion of parents (only matters with --widen):
        th["m"] = min(th["m"], 0.99)
    for x in args.fix:
        name, value = x.split("=")
        assert name in th, f"{name} isn't a parameter of {args.model}"
        th[name] = float(value)

    if crash:
        if args.ne_form == "linear":
            ne = th["floor"] + th["scale"] * trap
        else:
            ne = th["scale"] * (trap + 1) ** th["b"]
    else:
        ne = np.full(N_GEN, th["Ne_c"])
    ne = np.maximum(ne, 1.0)
    imm = "Ne_s" in th

    def mig(size):
        if "m" in th:
            return th["m"]
        return min(th["I"] / size, 0.99)

    # Generations (msprime time k) simulated with discrete-time
    # Wright-Fisher (see "ancestry model" below):
    small = [ne[N_GEN - 1 - k] < DTWF_MAX for k in range(N_GEN)]

    # ---- demography ----
    dem = msprime.Demography()
    # initially_active: SN is also the ancestral population of the split.
    dem.add_population(name="SN", initial_size=ne[N_GEN - 1],
                       initially_active=True)
    # In Wright-Fisher generations the source's size is capped at
    # GHOST_DTWF_MAX: Wright-Fisher's cost grows with population size (the
    # source can reach 1.1e8), and at >= 1e5 lineages that moved to the
    # source coalesce there with probability < 78 / 2e5 per pair over the
    # whole trap record, as under the uncapped size.
    def ghost_size(k):
        if small[k]:
            return min(th["Ne_s"], GHOST_DTWF_MAX)
        return th["Ne_s"]
    if imm:
        dem.add_population(name="ghost", initial_size=ghost_size(0))
        dem.set_migration_rate(source="SN", dest="ghost",
                               rate=mig(ne[N_GEN - 1]))
    for k in range(1, N_GEN):
        dem.add_population_parameters_change(
            time=k, initial_size=ne[N_GEN - 1 - k], population="SN")
        if imm and small[k] != small[k - 1]:
            dem.add_population_parameters_change(
                time=k, initial_size=ghost_size(k), population="ghost")
        if imm and "I" in th:
            dem.add_migration_rate_change(
                time=k, rate=mig(ne[N_GEN - 1 - k]), source="SN", dest="ghost")
    dem.add_population_parameters_change(
        time=N_GEN, initial_size=th["Ne_anc"], population="SN")
    if imm and small[N_GEN - 1]:
        dem.add_population_parameters_change(
            time=N_GEN, initial_size=th["Ne_s"], population="ghost")
    if imm and "I" in th:
        dem.add_migration_rate_change(
            time=N_GEN, rate=mig(th["Ne_anc"]), source="SN", dest="ghost")
    if args.older:
        dem.add_population_parameters_change(
            time=N_GEN + th["T_old"], initial_size=th["Ne_old"],
            population="SN")
    if imm:
        dem.add_population_split(time=N_GEN + th["T_s"], derived=["ghost"],
                                 ancestral="SN")
    dem.sort_events()

    samples = [msprime.SampleSet(n, population="SN", time=N_GEN - 1 - g,
                                 ploidy=2)
               for _, g, n in pools]

    # ---- ancestry model ----
    # Discrete-time Wright-Fisher (exact for small sizes) in generations
    # where SN's size is < DTWF_MAX, the standard coalescent elsewhere:
    # DTWF's cost grows with population size (5 s per locus at 1e7 vs.
    # 0.02 s for the coalescent), and the two agree when the size is much
    # larger than the sample (1,314 genomes).
    models = []
    k0 = 0
    for k in range(1, N_GEN + 1):
        if k == N_GEN or small[k] != small[k0]:
            if small[k0]:
                models.append(msprime.DiscreteTimeWrightFisher(duration=k - k0))
            else:
                models.append(msprime.StandardCoalescent(duration=k - k0))
            k0 = k
    models.append(msprime.StandardCoalescent())

    # ---- counts per pool, one locus at a time ----
    # (Loci simulated separately: in one sequence, unlinked loci made the
    # coalescent's cost grow with the square of the number of loci.)
    # Sample nodes are in SampleSet order, 2 per individual:
    n_samp = 2 * sum(n for _, _, n in pools)
    ind = np.zeros((n_samp, len(pools)), dtype=np.int32)
    i0 = 0
    for j, (_, _, n) in enumerate(pools):
        ind[i0:(i0 + 2 * n), j] = 1
        i0 += 2 * n

    W = args.locus_len
    all_pos = []
    all_counts = []
    for i in range(args.n_loci):
        ts = msprime.sim_ancestry(
            samples=samples, demography=dem, sequence_length=W,
            recombination_rate=args.rec, model=models,
            random_seed=int(rng.integers(1, 2**31 - 1)))
        ts = msprime.sim_mutations(
            ts, rate=args.mu, model=msprime.BinaryMutationModel(),
            discrete_genome=True, random_seed=int(rng.integers(1, 2**31 - 1)))
        if ts.num_sites == 0:
            continue
        G = ts.genotype_matrix().astype(np.int32)  # sites x samples
        all_counts.append((G @ ind).astype(np.int16))
        all_pos.append(i * W + ts.tables.sites.position.astype(np.int64))
    if all_counts:
        counts = np.concatenate(all_counts)
        out_pos = np.concatenate(all_pos).astype(np.int32)
    else:
        counts = np.zeros((0, len(pools)), dtype=np.int16)
        out_pos = np.zeros(0, dtype=np.int32)

    # Sites fixed for the derived allele in every sample aren't segregating:
    tot = np.array([2 * n for _, _, n in pools])
    seg = (counts.sum(axis=1) > 0) & ~np.all(counts == tot, axis=1)

    out_pos[seg].tofile(args.out + ".pos.bin")
    np.ascontiguousarray(counts[seg]).tofile(args.out + ".counts.bin")
    with open(args.out + ".params.tsv", "w") as f:
        f.write("param\tvalue\n")
        f.write(f"model\t{args.model}\n")
        f.write(f"seed\t{args.seed}\n")
        for x in sorted(th):
            f.write(f"{x}\t{th[x]:.8g}\n")
        f.write(f"n_sites\t{int(seg.sum())}\n")


if __name__ == "__main__":
    main()
