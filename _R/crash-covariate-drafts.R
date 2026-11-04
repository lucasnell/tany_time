#'
#' Draft crash covariates for BayPass, for the owner to choose from
#' (2026-11-04). Reads only abundances and the pool list: no allele
#' frequencies or BayPass results, so the choice can't depend on them.
#'
#' Generations are indexed in trap order (a then b each year). A crash
#' generation has trap N below CRASH_N; consecutive crash generations form
#' one crash episode. For each of the 18 time pools, the drafts are:
#'
#'   gens_since_crash  generations since the last generation of the most
#'                     recent crash episode (time to recover from it)
#'   recovering        1 if the pool was sampled after a crash episode and
#'                     before N first returned above RECOVERED_N, else 0
#'   min_logN_<k>      lowest logN in the k generations before the pool
#'                     (depth of the most recent bottleneck); k = 6, 10
#'   hm_logN_<k>       log of the harmonic mean of (N + 1) over the same k
#'                     generations (the long-run effective size of a
#'                     fluctuating population is near the harmonic mean)
#'
#' Also reported: each draft's correlation with logN across the 18 pools,
#' and the crash episodes for a range of CRASH_N values.
#'
#' Output: `_data/crash-covariate-drafts.csv`,
#' `_data/crash-covariate-drafts.txt`.
#'

library(tidyverse)

CRASH_N <- 100
RECOVERED_N <- 10000

abund <- read_csv("_data/tany-abundances.csv", show_col_types = FALSE) |>
    arrange(year, gen) |>
    mutate(g = row_number())

pools <- tibble(pool = read_lines("_bash/samples/time-pools.txt")) |>
    mutate(gen = if_else(str_sub(pool, 4, 4) == "1", "a", "b"),
           yy = as.integer(str_sub(pool, 6, 7)),
           year = if_else(yy < 70, 2000L + yy, 1900L + yy)) |>
    left_join(abund, by = c("year", "gen"))
stopifnot(!anyNA(pools$g))

#' Crash episodes: runs of consecutive generations with N < crash_n.
episodes <- function(crash_n) {
    abund |>
        mutate(crash = N < crash_n,
               ep = cumsum(crash & !lag(crash, default = FALSE))) |>
        filter(crash) |>
        group_by(ep) |>
        summarize(start = min(g), end = max(g),
                  label = sprintf("%i%s-%i%s", year[1], gen[1],
                                  year[n()], gen[n()]),
                  .groups = "drop")
}
eps <- episodes(CRASH_N)

draft_df <- pools |>
    rowwise() |>
    mutate(
        last_end = max(eps$end[eps$end < g], -Inf),
        gens_since_crash = g - last_end,
        recovered_since = any(abund$N[abund$g > last_end & abund$g < g] >
                                  RECOVERED_N),
        recovering = as.integer(!recovered_since & N <= RECOVERED_N),
        min_logN_6 = min(abund$logN[abund$g >= g - 6 & abund$g < g]),
        min_logN_10 = min(abund$logN[abund$g >= g - 10 & abund$g < g]),
        hm_logN_6 = log(1 / mean(1 / (abund$N[abund$g >= g - 6 &
                                                  abund$g < g] + 1))),
        hm_logN_10 = log(1 / mean(1 / (abund$N[abund$g >= g - 10 &
                                                   abund$g < g] + 1)))) |>
    ungroup() |>
    select(pool, year, gen, N, logN, gens_since_crash, recovering,
           starts_with("min_logN"), starts_with("hm_logN"))

write_csv(draft_df, "_data/crash-covariate-drafts.csv")

sink("_data/crash-covariate-drafts.txt", split = TRUE)
cat(sprintf("Crash episodes (N < %i):\n", CRASH_N))
print(select(eps, label), n = Inf)
cat("\nCrash episodes by threshold:\n")
for (cn in c(25, 50, 100, 150, 200, 250, 500)) {
    cat(sprintf("  N < %3i: %s\n", cn,
                paste(episodes(cn)$label, collapse = ", ")))
}
cat("\nDrafts per pool:\n")
print(mutate(draft_df, across(where(is.double) & !c(year, N), \(x) round(x, 2))),
      n = Inf, width = Inf)
cat("\nCorrelation with logN across the 18 pools:\n")
draft_df |>
    summarize(across(gens_since_crash:hm_logN_10, \(x) cor(x, logN))) |>
    pivot_longer(everything(), names_to = "draft", values_to = "r_logN") |>
    mutate(r_logN = round(r_logN, 2)) |>
    print()
sink()
