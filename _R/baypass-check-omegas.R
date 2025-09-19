
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
    discard(\(x) str_detect(x, "-sub0")) |>
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

# # A tibble: 3 × 3
#       a     b   fmd
#   <int> <int> <dbl>
# 1     3     1 3.78
# 2     3     2 3.77
# 3     2     1 0.272


#'
#' Because sub-sample 3 is quite dissimilar to the other ones (FMD > 3.7),
#' then we have to use the single covariance matrix
#'


#'
#' I'm going to use the one based on all SNPs, but out of curiosity,
#' let's check for how similar these sub-samples are
#' to the omega matrix calculated using all SNPs.
#'
omega0 <- paste0(baypass_init_dir,
                 "/baypass-init-sub0/tany-init-sub0_mat_omega.out") |>
    read_table(col_names = paste0("samp_", 1:17), col_types = cols()) |>
    as.matrix() |>
    unname()

map_dbl(1:3, \(i) fmd.dist(omega0, omegas[[i]]))
# [1] 0.5697506 0.6034402 4.1120763

#'
#' Looks like it's #3 that's the problem. < shrug >
#'
