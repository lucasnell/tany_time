# ABC design: crashes, immigration, and neutral diversity

Status: approved by the owner, 2026-11-12, before any observed summary
statistic was computed or any descriptive diversity result was read (TODO
review step 5). Changes after approval are listed at the end with dates
and reasons.

Question: are changes in neutral diversity across the SN time series
explained by local size changes that follow the trap counts, or do they
also need immigration from an unsampled source?

Method: ABC by random forests (`abcrf`; Pudlo et al. 2016,
doi:10.1093/bioinformatics/btv684; Raynal et al. 2019,
doi:10.1093/bioinformatics/bty867). Simulations use msprime (Hudson
coalescent; Baumdicker et al. 2022, doi:10.1093/genetics/iyab229).
Summary statistics for observed and simulated data come from one R
function.

## Time and data

- Generation g = 2 × (year − 1977) + (0 for `a`, 1 for `b`), so g runs
  from 0 (1977a) to 77 (2015b). Trap counts N_g are from
  `_data/tany-abundances.csv`.
- Samples: the 18 pools in `_bash/samples/time-pools.txt`. Each takes
  2 × n_adults genomes (from `_data/full-sample-info.csv`) at its
  generation.
- Mutation rate μ = 2.1e-9 per site per generation, fixed: direct
  estimate for *C. riparius* from mutation-accumulation lines (95% CI
  1.4e-9 to 3.1e-9; Oppold & Pfenninger 2017, Evolution Letters 1:86,
  doi:10.1002/evl3.8; Zotero `IARTHK3W`). The *D. melanogaster* rate
  (5.49e-9; Schrider et al. 2013) is a sensitivity run.
- Recombination rate r = 1e-8 per site per generation within loci, fixed
  (owner, 2026-11-12: rates fixed rather than free, to keep the model
  simple). No per-generation estimate exists for a chironomid. The *C. riparius*
  estimates are population-scaled (ρ = 4 Ne r) per bp:
  - 0.04–0.07 in 50-kb windows (LDhelmet; Schmidt et al. 2020, G3
    10:1151, doi:10.1534/g3.119.400710), "within the range of those
    estimated for *Drosophila melanogaster* (0.01 to 0.11)";
  - chromosome means 0.009–0.016 in 10-kb windows (iSMC; Pettrich &
    Waldvogel 2025, BMC Genomics 26:1002,
    doi:10.1186/s12864-025-12130-7, Table 1).

  Converting either needs *C. riparius* Ne, which neither paper gives.
  For ρ of 0.009 to 0.07 and Ne of 1e5 to 1e7, r runs from 2e-10 to 2e-7;
  *D. melanogaster* (similar ρ) is about 1e-8. The summary statistics are
  means over windows, which recombination within loci doesn't change, so r
  mainly sets the variance among windows.
- θ = 0.0089 (`00-references/04-ont-pi`), so θ / (4μ) ≈ 1.06e6. This is
  the reference value for the ancestral and source size priors.

## Models

Local effective size in generation g:

- crash-driven: Ne_g = floor + scale × N_g;
- constant: Ne_g = Ne_c.

Immigration comes from an unsampled source population ("ghost") of
constant size Ne_s. It split from SN T_s generations before 1977. Going
forward in time, a proportion m_g of SN's parents in generation g are
immigrants. In models with immigration, it also runs before 1977, at m
(M2, M3) or I / Ne_anc (M5), so SN and the source are near
migration–drift balance in 1977. The owner reports that four nearby lakes
are well connected to Mývatn (F_ST close to 0) and at times have large
populations, though far smaller than Mývatn's.

| ID | Local size | Immigration |
|---|---|---|
| M1 | crash-driven | none (m = 0) |
| M2 | crash-driven | constant proportion: m_g = m |
| M3 | constant | constant proportion: m_g = m |
| M4 | constant | none |
| M5 | crash-driven | constant number of immigrants per generation, I: m_g = min(I / Ne_g, 0.99), so the immigrant share rises during crashes |

- Before 1977, SN has constant size Ne_anc.
- Census size: the owner's rough estimate for Mývatn reaches 1.3e13
  adults. With a peak trap count of 181,832 (2000a), that is about 7e7
  adults per trap individual at the peak. Effective size per generation
  is bounded by Ne / N ≤ 0.14 at the peak (scale ≤ 1e7).
- Main contrast: M2 + M5 vs M1.
- Immigration as a group: M2 + M3 + M5 vs M1 + M4.

## Priors

All log-uniform and independent.

| Parameter | Models | Range | Reasoning |
|---|---|---|---|
| scale (Ne per trap individual) | M1, M2, M5 | 1 to 1e7 | at the peak trap count (181,832), Ne from 1.8e5 to 1.8e12, i.e. Ne / N from ~1e-8 to 0.14 of the census estimate (1.3e13) |
| floor (Ne when N = 0) | M1, M2, M5 | 10 to 1e9 | a zero trap count at SN needn't mean few adults lake-wide; from a near-total crash to crashes that barely reduce Ne |
| Ne_c | M3, M4 | 1e3 to 1e12 | brackets the crash-driven range |
| Ne_anc | all | 1.1e5 to 1.1e7 | θ / (4μ) ÷ 10 to × 10; not raised with the census, since observed π ≈ 4 Ne_anc μ pins it |
| Ne_s | M2, M3, M5 | 1.1e5 to 1.1e8 | as Ne_anc, up to 100×: nearby lakes may peak while Mývatn crashes (owner) |
| T_s (generations) | M2, M3, M5 | 1e2 to 1e7 | with immigration before 1977, the split time matters only when m is small; kept wide |
| m | M2, M3 | 1e-5 to 0.5 | from negligible to half of parents |
| I (immigrants per generation) | M5 | 1 to 1e9 | scaled with the local Ne range; the share I / Ne_g is capped at 0.99 |

## Summary statistics

**Data.** 2,000 random 10-kb windows from the partially masked syncs.
Each window must have ≥ 5,000 sites unmasked in all 18 pools. Simulated
locus i takes the observed per-pool depth and masking of window i, site by
site.

**Reads in simulations.**
- Reference reads at each site ~ Binomial(depth, p), where p is the
  sample frequency in the pool.
- Errors: each read is flipped with probability e = 0.00023, estimated
  from reads of a third or fourth base at the observed sites (pools
  0.00019-0.00029; `_data/abc-inputs.txt`); was 0.002 (see "Changes after
  approval").
- Filters as in `01-snape` masking: coverage ≥ 10 and ≤ the pool's 95th
  percentile, already applied through the observed masks.
- A site counts as segregating in a pool if the minor allele has ≥ 2
  reads.

**Per pool** (18 pools × 4 = 72 statistics), as means over windows:
1. π, corrected for pool size n and depth c:
   (c / (c − 1)) × (n / (n − 1)) × 2p(1 − p), averaged over sites.
2. Watterson's θ with the pool-seq correction of Kofler et al. 2011
   (PoPoolation; doi:10.1371/journal.pone.0015925).
3. Tajima's D from items 1 and 2, as in Kofler et al. 2011 (now
   π / θ_W − 1; see "Changes after approval").
4. Rare-allele share: among segregating sites, the proportion with minor
   allele frequency < 0.05.

**Temporal** (34 statistics): F between consecutive pools (17) and from
SN-2-92 to each later pool (17). Uses the Jónás et al. 2016 sampling
correction (plan II; doi:10.1534/genetics.116.191197), on sites
segregating in either pool.

**Among pools** (153 statistics): pairwise F_ST with the pool-seq
estimator of Hivert et al. 2018 (doi:10.1534/genetics.118.300900; as in
`poolfstat`).

`abcrf` adds linear discriminant axes of these 259 statistics by default.

## Simulations

- 25,000 per model (125,000 in total) for the reference table, plus 1,000
  per model held out for checking parameter recovery.
- Pilot first: 100 per model, to measure run time and run the health
  checks in the plan. If 25,000 per model isn't feasible, use 10,000 per
  model at least.
- Posterior-predictive check: 1,000 simulations from the posterior of the
  most probable model.

## Decision rule

- Report the posterior probability of every model.
- Report the immigration group (M2 + M3 + M5) vs the no-immigration group
  (M1 + M4), with the out-of-bag prior error rate for both the
  five-model and the two-group choice.
- Immigration is called "supported" only if the immigration group's
  posterior is ≥ 0.9.
- Parameters (scale, floor, m or I, Ne_s, T_s) are estimated with
  `regAbcrf` under each model with posterior ≥ 0.1, reporting medians and
  90% intervals.
- m or I is called identifiable only if 90% intervals cover the true
  values in 85–95% of held-out simulations and are narrower than the
  prior.

## Sensitivity (all reported, none selected)

1. Local size as Ne_g = scale × (N_g + 1)^b instead of floor + scale ×
   N_g, with b log-uniform on 0.1 to 1.
2. All prior ranges widened 10× on each side.
3. Leaving out SN-1-95 (included by owner override), SN-1-13 (n = 14), and
   SN-2-97 (n = 19).
4. Per-pool diversity statistics alone vs temporal and among-pool
   statistics alone.
5. Older history: one extra size change before 1977, with a size
   log-uniform from 1.1e4 to 1.1e8 and a time log-uniform from 1e2 to 1e5
   generations before 1977.
6. Mutation rate 5.49e-9 (*D. melanogaster*), with Ne_anc and Ne_s
   priors re-centred on θ / (4μ) ≈ 4.05e5 by the same rules (4e4 to 4e6;
   4e4 to 4e7).

## Changes after approval

2026-11-12, while writing the code (`_R/abc-inputs.R`,
`_bash/bin/abc-sim.py`, `_bash/bin/abc-stats.R`); none uses observed
statistics, which haven't been computed. Items 1-2 and 9 change approved
choices (approved by the owner, 2026-11-12); 3-8 are implementation
details.

1. Tajima's D is replaced by π / θ_W − 1 (same sign; no variance
   normalization). The pool-seq variance of D (Kofler et al. 2011) isn't
   implemented; since the same code runs on observed and simulated data,
   the normalization adds no information for the ABC.
2. Temporal F uses sites covered in both pools and segregating (minor
   allele ≥ 2 reads) in either, not Jónás et al.'s rule (polymorphic in the
   first sample). With no drift, their rule gave F of −0.012 to −0.0007
   (bias from conditioning on the first sample); this rule gives +0.0001
   to +0.0038, largest at low depth (`_bash/tests/abc-stats-check.R` →
   `_data/abc-stats-check.txt`). Either bias cancels in the ABC; it
   matters for F converted to Ne.
3. 258 statistics, not 259: the first-to-second pool pair is already a
   consecutive pair, so there are 16 first-to-later pairs.
4. Ancestry model: discrete-time Wright-Fisher (exact for small sizes) in
   generations where SN's size is < 1e4, the standard coalescent
   otherwise. Wright-Fisher's cost grows with population size (5 s per
   locus at 1e7 vs. 0.02 s for the coalescent); the two agree when the
   size is much larger than the sample (1,314 genomes). In Wright-Fisher
   generations the source's size is capped at 1e5 (lineages there coalesce
   with probability < 78 / 2e5 per pair over the trap record either way).
5. Loci are simulated one at a time (in one sequence, unlinked loci made
   the coalescent's cost grow with the square of the number of loci).
6. Sites that don't segregate but show ≥ 2 error reads of one base are
   simulated too, as the observed data contain them.
7. Depths > 255 are treated as masked in both observed and simulated
   data (depths are stored as one byte).
8. Observed alleles at a site: the two bases with the most reads summed
   over pools; reads of the other two bases are dropped and used to
   estimate the sequencing error rate (`_data/abc-inputs.txt`).
9. Sequencing error rate 0.00023 instead of 0.002: the rate estimated
   from third- and fourth-base reads in the observed windows
   (`_R/abc-inputs.R` → `_data/abc-inputs.txt`; pools 0.00019-0.00029).
   With 0.002, simulations would carry ~9x more error reads than the
   data, inflating θ_W and the rare-allele share in simulations only.
