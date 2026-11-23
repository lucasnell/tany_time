#'
#' Which observed statistics do the models without immigration (M1, M4)
#' fail to reach? (2026-11-23; diagnostic of the main ABC result in
#' `_data/abc-results.txt`, not a new test.) The main run gives
#' P(immigration) = 1, while 27-30 of the 235 statistics are outside the
#' central 98% of M1's and M4's simulations and inside those of M2, M3,
#' and M5. If those few statistics carry the model choice, the result
#' rests on a misfit of the no-immigration models in those statistics
#' rather than on the data as a whole.
#'
#' Per statistic, the observed value's quantile among each model's
#' simulations (`_R/abc-analysis.R` reference table, variant `main`),
#' and for the statistics outside M1's and M4's central 98%, their type
#' (pi, tw, tdr, rare, fc, fst), pools, and gap in generations.
#'
#' Writes `_data/abc-misfit.txt`.
#'

library(tidyverse)

source("_R/data-dir.R")

out_file <- "_data/abc-misfit.txt"
sink(out_file, split = TRUE)

abc_dir <- tany_dir("abc")
obs_df <- read_tsv(file.path(abc_dir, "abc-observed-stats.tsv"),
                   show_col_types = FALSE)
obs <- setNames(obs_df$value, make.names(obs_df$stat))
obs <- obs[!str_detect(names(obs), fixed("SN.1.95_S59"))]
stopifnot(length(obs) == 235)
st <- names(obs)

ref <- list.files(file.path(abc_dir, "sims"), "^abc-main-M[1-5]-b\\d+\\.tsv\\.gz$",
                  full.names = TRUE) |>
    map(\(f) read_tsv(f, show_col_types = FALSE, progress = FALSE,
                      col_types = cols(model = "c", .default = "d"))) |>
    list_rbind()
names(ref) <- make.names(names(ref))
stopifnot(nrow(ref) == 50000)

q <- map(paste0("M", 1:5), \(md) {
    x <- ref[ref$model == md, st]
    tibble(stat = st, model = md, q = map2_dbl(x, obs, \(s, o) mean(s <= o)))
}) |>
    list_rbind() |>
    pivot_wider(names_from = model, values_from = q)

gen_of <- function(p) {
    y <- as.integer(str_sub(p, 6, 7))
    y <- ifelse(y < 70, 2000 + y, 1900 + y)
    2 * (y - 1977) + as.integer(str_sub(p, 4, 4)) - 1
}
pp <- str_match(q$stat, "^(fc|fst)_(SN\\.\\d\\.\\d{2}_S\\d+)_(SN\\.\\d\\.\\d{2}_S\\d+)$")
q <- q |>
    mutate(type = str_extract(stat, "^[a-z]+"),
           gap = abs(gen_of(pp[, 3]) - gen_of(pp[, 4])),
           out_none = pmin(M1, M4) > 0.99 | pmax(M1, M4) < 0.01,
           in_imm = pmin(M2, M3, M5) >= 0.01 & pmax(M2, M3, M5) <= 0.99)

cat("== Statistics with the observed value outside the central 98% of both",
    "M1 and M4 (on the same side)\n")
cat("Count:", sum(q$out_none), "of", nrow(q), "; of these, inside the",
    "central 98% of all of M2, M3, and M5:", sum(q$out_none & q$in_imm), "\n")
cat("\nBy type (all statistics vs. outside M1 and M4):\n")
q |>
    group_by(type) |>
    summarize(n = n(), outside_none = sum(out_none), .groups = "drop") |>
    as.data.frame() |>
    print(row.names = FALSE)
cat("\nThose statistics (observed quantile among each model's simulations;",
    "side = whether the observed value is above or below M1 and M4):\n")
q |>
    filter(out_none) |>
    mutate(side = ifelse(M1 > 0.99, "above", "below"),
           across(M1:M5, \(x) round(x, 3))) |>
    arrange(type, gap) |>
    select(stat, type, gap, side, M1:M5) |>
    as.data.frame() |>
    print(row.names = FALSE)

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
