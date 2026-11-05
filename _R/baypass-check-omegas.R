#'
#' Check consistency of Omega (population covariance) matrices from
#' `01-baypass-init`: 3 thinned subsamples for each of 3 missing-data
#' settings (2026-11 rerun).
#'
#' Within a setting, the cluster step `02-check-omegas` already reports
#' pairwise FMD; this recomputes it (should match `omega-check.tsv`) and
#' adds pairs across settings. Settings differ in their SNP sets, so FMD
#' across settings shows how much the missing-data setting itself changes
#' Omega, compared with the subsample-to-subsample noise within a setting.
#'
#' Reads the `baypass-init-<tag>-sub<k>.tar.gz` outputs in the local data
#' folder (`tany_dir("baypass-init")`). Writes `_data/baypass-omega-fmd.csv`.
#'

library(tidyverse)
library(geigen) # geigen function inside fmd.dist

source("_R/data-dir.R")


init_dir <- tany_dir("baypass-init")
mp_tags <- c("mp050", "mp075", "mp100")



# =============================================================================*
# ---- read Omegas ----
# =============================================================================*


omega_df <- crossing(tag = mp_tags, sub = 1:3) |>
    mutate(name = sprintf("%s-sub%i", tag, sub))

#' Extract only the Omega file from each archive.
omegas <- omega_df$name |>
    set_names() |>
    map(\(nm) {
        fn <- sprintf("baypass-init-%s/tany-init-%s_mat_omega.out", nm, nm)
        untar(file.path(init_dir, sprintf("baypass-init-%s.tar.gz", nm)),
              files = fn, exdir = tempdir())
        file.path(tempdir(), fn) |>
            read_table(col_names = FALSE, col_types = cols()) |>
            as.matrix() |>
            unname()
    })

# Should be TRUE (all 18 x 18 and symmetric):
all(map_lgl(omegas, \(x) isSymmetric(x) && all(dim(x) == 18)))



# =============================================================================*
# ---- FMD ----
# =============================================================================*


#' From baypass `utils` folder
#'
#' This function computes the metric proposed by Forstner and Moonen (2003)
#' to evaluate the distance between two covariance matrices (FMD distance).
#'
#' From Gautier 2015 (doi: 10.1534/genetics.115.181453),
#' ideally all pairwise FMD are < 0.5
#'
fmd.dist <- function(mat1, mat2) {
    return(sqrt(sum(log(geigen(mat1,mat2)$values)**2)))
}

fmd_df <- combn(omega_df$name, 2) |>
    t() |>
    as_tibble(.name_repair = \(x) c("a", "b")) |>
    mutate(tag_a = str_remove(a, "-sub.*"),
           tag_b = str_remove(b, "-sub.*"),
           pair = ifelse(tag_a == tag_b, "within", "across"),
           fmd = map2_dbl(a, b, \(i, j) fmd.dist(omegas[[i]], omegas[[j]])))

write_csv(fmd_df, "_data/baypass-omega-fmd.csv")


#' Within-setting pairs should match the cluster's `omega-check.tsv`.
fmd_df |>
    filter(pair == "within") |>
    select(a, b, fmd)

fmd_df |>
    group_by(pair, tag_a, tag_b) |>
    summarize(min = min(fmd), median = median(fmd), max = max(fmd),
              .groups = "drop")

#   pair   tag_a tag_b   min median   max
# 1 across mp050 mp075 0.283  0.320 0.391
# 2 across mp050 mp100 0.557  0.593 0.637
# 3 across mp075 mp100 0.432  0.456 0.518
# 4 within mp050 mp050 0.222  0.282 0.288
# 5 within mp075 mp075 0.254  0.288 0.318
# 6 within mp100 mp100 0.285  0.330 0.334

#'
#' (2026-11-05) Within-setting values match `omega-check.tsv`. mp050 and
#' mp075 Omegas differ about as much as subsamples within a setting do;
#' mp100 (data required in all pools) differs from mp050 by more
#' (0.56-0.64) than any within-setting pair, so the missing-data setting
#' changes Omega, not just the SNP count.
#'
