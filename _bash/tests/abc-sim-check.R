#'
#' Health checks of the ABC simulator and statistics together
#' (`_bash/bin/abc-sim.py` → `abc-stats.R`), run locally from the repo
#' root with a python that has msprime 1.4.4 and tskit 1.0.3:
#'   ABC_PYTHON=<python> Rscript _bash/tests/abc-sim-check.R
#' Output: `_data/abc-sim-check.txt`.
#'
#' Real pools (sizes and generations) and trap counts; stand-in depths
#' (Poisson, mean 60, 20% of sites masked); 200 loci; no sequencing
#' error, except in check 4. Parameters fixed with `--fix`.
#'   1. Constant size 1e6 throughout (M4, Ne_c = Ne_anc): π and θ_W ≈
#'      4 Ne μ = 0.0084; pi / thw - 1 ≈ 0; temporal F ≈ 0 (+ the small
#'      depth bias in `_data/abc-stats-check.txt`).
#'   2. Size 1e3 from 1977 (M4, Ne_c = 1e3, Ne_anc = 1e6): temporal F
#'      over t generations ≈ 1 - exp(-t / 2000).
#'   3. Immigration from a diverged source (M3, Ne_c = 1e3, m = 0.1, T_s
#'      = 1e6, Ne_s = 1e6): F stays far below check 2, and F_ST among pools
#'      would stay near 0 (F_ST needs poolfstat, NA locally).
#'   4. As check 1 with sequencing error e = 0.002.
#'

source("_bash/bin/abc-stats.R")
Sys.setenv(ABC_SIM_PY = "_bash/bin/abc-sim.py")

out_file <- "_data/abc-sim-check.txt"
sink(out_file, split = TRUE)

L <- 200L
in_dir <- tempfile("abc-simcheck")
dir.create(in_dir)
set.seed(1263407051)

#' Pools and trap counts as in `_R/abc-inputs.R`:
suppressMessages(library(tidyverse))
pool_df <- read_csv("_data/full-sample-info.csv", show_col_types = FALSE) |>
    distinct(biotech_id, n_adults) |>
    filter(biotech_id %in% read_lines("_bash/samples/time-pools.txt")) |>
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
gens <- setup$pools$gen
pr <- setup$pairs_t
dt <- gens[pr[2, ]] - gens[pr[1, ]]

run_one <- function(model, fix, err = 0) {
    Sys.setenv(ABC_ERR = err)
    f <- tempfile(fileext = ".tsv")
    run_sims(in_dir, setup, model, seed0 = 1L, n_sims = 1L, out_file = f,
             extra = c(rbind("--fix", fix)))
    x <- read.delim(f)
    x
}

summ <- function(x, label) {
    s <- unlist(x[1, ])
    pi_ <- as.numeric(s[grepl("^pi_", names(s))])
    thw <- as.numeric(s[grepl("^thw_", names(s))])
    tdr <- as.numeric(s[grepl("^tdr_", names(s))])
    fc <- as.numeric(s[grepl("^fc_", names(s))])
    cat("\n== ", label, " (simulation ", s[["t_sim"]], " s, total ",
        s[["t_all"]], " s)\n", sep = "")
    cat("pi: mean ", signif(mean(pi_), 3), " (range ",
        signif(min(pi_), 3), "-", signif(max(pi_), 3), ")\n", sep = "")
    cat("thw: mean ", signif(mean(thw), 3), "; pi / thw - 1: mean ",
        signif(mean(tdr), 3), "\n", sep = "")
    invisible(fc)
}

fc1 <- summ(run_one("M4", c("Ne_c=1e6", "Ne_anc=1e6")),
            "1. Constant size 1e6 (expected pi = thw = 0.0084)")
cat("temporal F: mean ", signif(mean(fc1), 3), " (range ",
    signif(min(fc1), 3), " to ", signif(max(fc1), 3), ")\n", sep = "")

fc2 <- summ(run_one("M4", c("Ne_c=1e3", "Ne_anc=1e6")),
            "2. Size 1e3 from 1977")
fc3 <- summ(run_one("M3", c("Ne_c=1e3", "Ne_anc=1e6", "m=0.1", "T_s=1e6",
                            "Ne_s=1e6")),
            "3. Size 1e3 from 1977, immigration m = 0.1")
cat("\nTemporal F by pair (t generations; check 2 expected 1 - exp(-t / 2000)):\n")
print(data.frame(pair = names(fc1) <- paste(setup$pools$pool[pr[1, ]],
                                             setup$pools$pool[pr[2, ]],
                                             sep = " > "),
                 t = dt, check2_expected = signif(1 - exp(-dt / 2000), 2),
                 check2 = signif(fc2, 2), check3 = signif(fc3, 2)),
      row.names = FALSE)

summ(run_one("M4", c("Ne_c=1e6", "Ne_anc=1e6"), err = 0.002),
     "4. As check 1 with sequencing error 0.002")

sink()
writeLines(sub("[[:space:]]+$", "", readLines(out_file)), out_file)
