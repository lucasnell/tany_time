#'
#' Health check of the `older` sensitivity variant (`_R/abc-design.md`,
#' sensitivity item 5: one size change of SN before 1977), run locally
#' from the repo root before it's simulated on the cluster:
#'   ABC_PYTHON=<python> Rscript _bash/tests/abc-older-check.R
#' Output: `_data/abc-older-check.txt`.
#'
#' Setup as `abc-sim-check.R` (real pools and trap counts; stand-in depths;
#' 200 loci; sequencing error 0.00023 and the template step, with the same
#' seed in every run, so the same template draws). M4 with Ne_c = Ne_anc =
#' 1e6, and before that, from T_old generations before 1977, size Ne_old:
#'   1. no `--older` (reference): pi / thw - 1 near the M4 value;
#'   2. expansion: Ne_old = 1e5, T_old = 1e4 -> excess of rare alleles
#'      (rare share up, pi / thw - 1 down) and pi below check 1;
#'   3. contraction: Ne_old = 1e7, T_old = 1e4 -> the reverse;
#'   4. expansion at the prior's oldest time: Ne_old = 1e5, T_old = 1e5
#'      (still recent relative to 2 Ne_anc = 2e6 generations, so the
#'      excess should remain, though smaller than in check 2 or not).
#' Also checks that Ne_old and T_old are recorded as parameters.
#'

source("_bash/bin/abc-stats.R")
Sys.setenv(ABC_SIM_PY = "_bash/bin/abc-sim.py")
suppressMessages(library(tidyverse))

out_file <- "_data/abc-older-check.txt"
sink(out_file, split = TRUE)

L <- 200L
in_dir <- tempfile("abc-oldercheck")
dir.create(in_dir)
set.seed(1554902261)

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
file.copy("_bash/04-diversity/abc-priors.csv", in_dir)

setup <- abc_setup(in_dir)

run_one <- function(fix, older) {
    f <- tempfile(fileext = ".tsv")
    run_sims(in_dir, setup, "M4", seed0 = 1L, n_sims = 1L, out_file = f,
             extra = c(if (older) "--older", rbind("--fix", fix)))
    read.delim(f, check.names = FALSE)
}

base <- c("Ne_c=1e6", "Ne_anc=1e6")
runs <- list(
    "1. no older change" = run_one(base, FALSE),
    "2. expansion (Ne_old 1e5, T_old 1e4)" =
        run_one(c(base, "Ne_old=1e5", "T_old=1e4"), TRUE),
    "3. contraction (Ne_old 1e7, T_old 1e4)" =
        run_one(c(base, "Ne_old=1e7", "T_old=1e4"), TRUE),
    "4. expansion at the prior's oldest (Ne_old 1e5, T_old 1e5)" =
        run_one(c(base, "Ne_old=1e5", "T_old=1e5"), TRUE))

cat("Parameters recorded with --older:",
    all(c("Ne_old", "T_old") %in% names(runs[[2]])), "\n\n")
imap(runs, \(x, nm) {
    m <- \(p) mean(as.numeric(unlist(x[1, grepl(p, names(x))])))
    tibble(run = nm, pi = signif(m("^pi_"), 3), thw = signif(m("^thw_"), 3),
           tdr = round(m("^tdr_"), 3), rare = round(m("^rare_"), 3),
           t_all_s = x$t_all)
}) |>
    list_rbind() |>
    as.data.frame() |>
    print(row.names = FALSE)

sink()
writeLines(sub("[[:space:]]+$", "", readLines(out_file)), out_file)
