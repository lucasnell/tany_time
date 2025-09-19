
#'
#' Check consistency of Omega matrices across runs.
#'
#' First download, extract baypass-init runs into directory `baypass_init_dir`.
#'

library(tidyverse)
library(geigen) # geigen function inside fmd.dist


baypass_init_dir <- "~/_data/_baypass/init"


n_omegas <- list.dirs(baypass_init_dir) |>
    keep(\(x) str_detect(x, "-sub")) |>
    length()

omegas <- sprintf(paste0(baypass_init_dir, "/baypass-init-sub%i/",
                         "tany-init-sub%i_mat_omega.out"), 1:n_omegas, 1:n_omegas) |>
    map(\(fn) {
        fn |>
            read_table(col_names = paste0("samp_", 1:17), col_types = cols()) |>
            as.matrix() |>
            unname()
    })

# Should be true:
all(map_lgl(omegas, isSymmetric))


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

fmd_df <- crossing(a = 1:length(omegas), b = 1:length(omegas)) |>
    filter(a > b) |>
    mutate(fmd = map2_dbl(a, b, \(i, j) fmd.dist(omegas[[i]], omegas[[j]]))) |>
    arrange(desc(fmd))

range(fmd_df$fmd)
# [1] 0.3019383 0.3829098

fmd_df

# # A tibble: 10 × 3
#        a     b   fmd
#    <int> <int> <dbl>
#  1     5     2 0.383
#  2     3     2 0.367
#  3     2     1 0.356
#  4     5     3 0.340
#  5     4     2 0.339
#  6     5     1 0.329
#  7     4     3 0.325
#  8     3     1 0.321
#  9     4     1 0.320
# 10     5     4 0.302


#'
#' Because all are quite similar (all pairwise FMD < 0.4), then we randomly
#' choose an Omega matrix to use:
#'
set.seed(1991347408); sample.int(n_omegas, 1)
# [1] 3
