#'
#' Can the template-sampling step (`sim_reads()` in `_bash/bin/abc-stats.R`,
#' ABC design change 13) produce the differentiation between nearby pools
#' that failed the prior-predictive check (`_data/abc-results.txt`)? Run
#' locally from the repo root with a python that has msprime 1.4.4:
#'   ABC_PYTHON=<python> Rscript _bash/tests/abc-tmpl-check.R
#' Output: `_data/abc-tmpl-check.txt`.
#'
#' Setup as `_bash/tests/abc-sim-check.R` (real pools and trap counts,
#' stand-in depths with mean 60, 200 loci, the design's error rate), with
#' almost no drift (constant size 1e6, model M4), so any differentiation
#' between pools comes from sampling. The same number of templates per site
#' is given to every pool: none (the old model), 10, 30, 100, 300. Reported:
#' median temporal F and F_ST over pairs of pools 1-2 generations apart,
#' next to the observed range for those pairs (from `_data/abc-results.txt`
#' and the observed statistics); the same for F_ST (needs poolfstat).
#'

source("_bash/bin/abc-stats.R")
Sys.setenv(ABC_SIM_PY = "_bash/bin/abc-sim.py")
suppressMessages(library(tidyverse))

out_file <- "_data/abc-tmpl-check.txt"
sink(out_file, split = TRUE)

L <- 200L
in_dir <- tempfile("abc-tmplcheck")
dir.create(in_dir)
set.seed(1592744413)
pool_df <- read_csv("_data/full-sample-info.csv", show_col_types = FALSE) |>
    distinct(biotech_id, n_adults) |>
    filter(biotech_id %in% read_lines("_bash/samples/time-pools-18.txt")) |>
    mutate(yy = as.integer(str_match(biotech_id, "^SN-[12]-(\\d{2})")[, 2]),
           year = ifelse(yy < 70, 2000L, 1900L) + yy,
           g12 = as.integer(str_match(biotech_id, "^SN-([12])-")[, 2]),
           gen = 2L * (year - 1977L) + g12 - 1L) |>
    arrange(gen) |>
    select(pool = biotech_id, gen, n_adults)
write_tsv(pool_df, file.path(in_dir, "pools.tsv"))
read_csv("_data/tany-abundances.csv", show_col_types = FALSE) |>
    arrange(year, gen) |>
    getElement("N") |>
    write_lines(file.path(in_dir, "trap.txt"))
write_tsv(tibble(locus = 1:L, contig = "X", window = 1:L, start = 1),
          file.path(in_dir, "windows.tsv"))
P <- nrow(pool_df)
d <- rpois(W * P * L, 60)
d[runif(length(d)) < 0.2 | d < 10] <- 0L
writeBin(as.raw(pmin(d, MAX_DEP)), file.path(in_dir, "depth.bin"))
invisible(file.copy("_bash/04-diversity/abc-priors.csv", in_dir))

setup <- abc_setup(in_dir)
gens <- setup$pools$gen
pr <- setup$pairs_t
near <- abs(gens[pr[2, ]] - gens[pr[1, ]]) <= 2
fc_names <- paste("fc", setup$pools$pool[pr[1, ]], setup$pools$pool[pr[2, ]],
                  sep = "_")[near]
pf <- setup$pairs_f
near_f <- abs(gens[pf[2, ]] - gens[pf[1, ]]) <= 2
fst_names <- paste("fst", setup$pools$pool[pf[1, ]], setup$pools$pool[pf[2, ]],
                   sep = "_")[near_f]

obs <- read_tsv(file.path(Sys.getenv("ABC_OBS_DIR",
                                     "/Volumes/lan_images4/tany_time/abc"),
                          "abc-observed-stats.tsv"), show_col_types = FALSE)
obs_fc <- obs$value[obs$stat %in% fc_names]
obs_fst <- obs$value[obs$stat %in% fst_names]
cat("Observed F_ST, pairs <= 2 generations apart (", length(obs_fst),
    " pairs): median ", signif(median(obs_fst), 3), ", range ",
    signif(min(obs_fst), 3), "-", signif(max(obs_fst), 3), "\n", sep = "")
cat("Observed temporal F, pairs <= 2 generations apart (", length(obs_fc),
    " pairs, SN-1-95 pairs included): median ", signif(median(obs_fc), 3),
    ", range ", signif(min(obs_fc), 3), "-", signif(max(obs_fc), 3), "\n",
    sep = "")

res <- map(c(NA, 10, 30, 100, 300), \(tm) {
    setup$tmpl_prior <- if (is.na(tm)) NULL else c(tm, tm)
    f <- tempfile(fileext = ".tsv")
    Sys.setenv(ABC_ERR = 0.00023)
    run_sims(in_dir, setup, "M4", seed0 = 7L, n_sims = 1L, out_file = f,
             extra = c("--fix", "Ne_c=1e6", "--fix", "Ne_anc=1e6"))
    x <- read.delim(f, check.names = FALSE)
    v <- unlist(x[1, fc_names])
    w <- unlist(x[1, fst_names])
    tibble(templates = ifelse(is.na(tm), "none", as.character(tm)),
           median_fc = median(v), min_fc = min(v), max_fc = max(v),
           median_fst = median(w), min_fst = min(w), max_fst = max(w))
}) |>
    list_rbind()
cat("\nSimulated temporal F and F_ST, same pairs (constant size 1e6, no",
    "immigration):\n")
print(as.data.frame(res |> mutate(across(where(is.double), \(x) signif(x, 3)))))

sink()
writeLines(sub("[[:space:]]+$", "", readLines(out_file)), out_file)
