# TODO

## Housekeeping

- [ ] `_R/make-time-series-plot.R` saves its main figure to
      `~/Desktop/freq-abund-time.pdf`; point it at `_plots/`.
- [ ] Decide whether to move legacy root-level scripts into a subfolder.
- [ ] Move BLAST inputs/outputs off `~/Desktop` (`Tgraci_proteins.faa`,
      `Tgraci-funct-Dmel.csv`) to a stable path.
- [ ] `baypass-check-omegas.R`: the comment `range(fmd_df$fmd)` → 0.30–0.38
      contradicts the printed table (max 3.78). Update the stale comment.
- [ ] `README.md` contains only "midge_archive".

## Analysis questions

- [ ] Ising β. Moved from β = 1 (commit 6df3656) to β = 0 (commit 0ef54a3)
      because the assembly is not chromosome-level. Note that the Ising
      prior only links SNPs adjacent in the input file, so unknown linkage
      between contigs is ignored under either β; the β = 1 risk specific to
      this assembly is the ~44 contig junctions being treated as neighbours.
      Consider running β = 1 as a sensitivity check and reporting overlap
      with the β = 0 candidates. Only `beta0` outputs are local.
- [ ] Candidate stability: every candidate passes 30 dB in only 5–7 of 10
      runs. Report this alongside the median-based list rather than treating
      the 33 as a fixed set.
- [ ] Bayes factors hit a ceiling (52.96 dB = PIP of 1.0 across all 1,000
      MCMC samples). Top candidates are tied, so ordering among them by
      `bf_db` carries no information. More samples (`-nval`) would raise
      the ceiling.
- [ ] GO enrichment (`enricher`) uses the default universe — every
      *D. melanogaster* protein with GO annotation — not the genes that could
      have been hit (Tgraci genes near tested SNPs with a Dmel ortholog). With
      8 input genes, each significant term is probably driven by one gene.
      Needs a matched universe before any enrichment result is used.
- [ ] Allele polarity: BayPass uses reference-allele counts; `timeseries-data.csv`
      uses SNAPE minor-allele frequency. Check sign before interpreting
      direction of association.
- [ ] Time-series analysis for balancing selection on
      `timeseries-data.csv` — not started in this repo.
