
library(Rcpp)
library(tidyverse)



snape_file <- function(x){
    paste0("~/_data/_tany-time/snape/snape_masked_time/", x)
}
count_file <- function(x){
    paste0("~/_data/_tany-time/snape/count_masked_time/", x)
}


Rcpp::sourceCpp("_R/snape-helpers.cpp")



# =============================================================================*
# Read data ----
# =============================================================================*

samp_df <- read_csv("_data/full-sample-info.csv", col_types = "cDccidd")

pools <- read_lines(snape_file("snape_masked_time.names.gz"))
n_adults <- pools |>
    set_names() |>
    map_int(\(x) samp_df$n_adults[samp_df$biotech_id == x]) |>
    as.list()

n_alleles <- read_tsv(count_file("count_masked_time_noblanks.sync.gz"),
         col_names = c("contig", "pos", "ref", pools),
         col_types = paste(c("cic", rep("c", length(pools))),
                           collapse = ""),
         progress = FALSE) |>
    mutate(across(all_of(pools), split_count_strings)) |>
    select(all_of(pools)) |>
    as.list() |>
    count_alleles()

snape_df <- read_tsv(snape_file("snape_masked_time_noblanks.sync.gz"),
                     col_names = c("contig", "pos", "ref", pools),
                     col_types = paste(c("cic", rep("c", length(pools))),
                                       collapse = ""),
                     progress = FALSE) |>
    mutate(contig = factor(contig, levels = paste0("CONTIG", 1:45)),
           across(all_of(pools), get_snape_af)) |>
    filter(n_alleles == 2) |>
    # Kapun et al. (2021) only used SNPs > 500 bp away from each other.
    filter(min_dist_filter(contig, pos, 500)) |>
    # Filtering out softmasked sequences bc SNAPE script doesn't handle these
    # properly (yet)
    filter(!ref %in% letters)
count_df <- read_tsv(count_file("count_masked_time_noblanks.sync.gz"),
                     col_names = c("contig", "pos", "ref", pools),
                     col_types = paste(c("cic", rep("c", length(pools))),
                                       collapse = ""),
                     progress = FALSE) |>
    mutate(contig = factor(contig, levels = paste0("CONTIG", 1:45)),
           across(all_of(pools), split_count_strings)) |>
    filter(n_alleles == 2) |>
    # Kapun et al. (2021) only used SNPs > 500 bp away from each other.
    filter(min_dist_filter(contig, pos, 500)) |>
    # Filtering out softmasked sequences bc SNAPE script doesn't handle
    # these properly (yet)
    filter(!ref %in% letters)




# =============================================================================*
# Combine datasets ----
# =============================================================================*


#' Below, `n_eff` is the "effective number of observations at a given locus,
#' conditional on read depth and number of chromosomes in the sample"
#' (Feder et al. 2012; doi: 10.1371/journal.pone.0048588)

freq_df <- snape_df |>
    pivot_longer(all_of(pools), names_to = "pool", values_to = "freq") |>
    full_join(count_df |>
                  mutate(across(all_of(pools), \(x) map_int(x, sum))) |>
                  pivot_longer(all_of(pools), names_to = "pool",
                               values_to = "n_reads"),
              by = c("contig", "pos", "ref", "pool")) |>
    mutate(n_haps = 2L * map_int(pool, \(p) n_adults[[p]]),
           n_eff = (n_reads * n_haps - 1) / (n_reads + n_haps))


#' Comparison of different types of files and compression (with larger dataset):
#'
#' file  comp    write       read        size
#' csv   none    3.3sec      2sec        430MB
#' csv   gz      9.2sec      3.4sec      64.5MB
#' csv   bz2     25.6sec     9.3sec      39.3MB
#' csv   xz      4min        5sec        31.4MB
#' rds   bz2     1.1min      11.2sec     28.5MB
#' rds   gz      1.4min      4.6sec      27.1MB
#' rds   xz      1.4min      4.7sec      27.0MB
#'
#' Based on this, I'll save my outputs as `rds` files with `gz` compression.


write_rds(freq_df, "_data/allele-freqs-time.rds", compress = "gz")




# =============================================================================*
# Side note ----
# =============================================================================*
#'
#' Below shows that there can seemingly be multiple different alternative
#' alleles, but this is not actually the case (since we've already filtered
#' for this with the `n_alleles` object).
#' What's happening is that my `which_allele` function sometimes thinks
#' the reference is the alternative when the minor-allele frequency estimate
#' is near 0.5.
#' So don't bother using `which_allele` and safely assume these frequencies
#' refer to the same minor allele type (since it never refers to the reference).
#'
#' ----------------------------------------------------------------------------*
#
# # Which allele does each frequency refer to?
# allele_i <- which_allele(count_df |> select(all_of(pools)) |> as.list(),
#                          snape_df |> select(all_of(pools)) |> as.matrix())
# unq_alleles <- apply(allele_i, 1, \(x) length(unique(x)))
# unq_alleles |> table()
# # unq_alleles
# #      1      2
# # 321366 110957
#
# # Now we figure out why there are two alternative minor alleles.
#
# # This tells us that the most common minor allele is never the
# # reference allele.
# is_it_ref <- map_lgl(1:nrow(count_df),
#                      \(i) {
#                          ref_i <- match(count_df[["ref"]][[i]], c("A","T","C","G"))
#                          minor_i <- table(allele_i[i,]) |>
#                              sort(decreasing=TRUE) |>
#                              base::`[`(1) |>
#                              names() |> as.integer()
#                          return(ref_i == minor_i)
#                      })
# sum(is_it_ref)
# # [1] 0
#
#
# ea_counts <- count_df |> filter(unq_alleles == 2)
# ea_snapes <- snape_df |> filter(unq_alleles == 2)
#
# get_ea_info <- function(ea_counts, ea_snapes) {
#     t0 <- Sys.time()
#     one_ea_info <- function(i) {
#         a <- ea_counts[i,pools] |> as.list() |> do.call(what = c) |>
#             do.call(what = rbind) |>
#             base::`[`(,1:4) |>
#             (\(x) { x / (rowSums(x) |> matrix(nrow(x), ncol(x))) })()
#         b <- ea_snapes[i,pools] |> unlist()
#         wa <- max.col(-abs(a - matrix(b, nrow(a), ncol(a))), "first")
#         st <- sort(table(wa))[1]
#         st_i <- st |> names() |> as.integer()
#         ref_i <- match(ea_counts[["ref"]][[i]], c("A","T","C","G"))
#         out <- tibble(is_ref = st_i == ref_i,
#                       freq_min = min(b[wa == st_i]),
#                       freq_max = max(b[wa == st_i]))
#         return(out)
#     }
#     out <- lapply(1:nrow(ea_counts), one_ea_info) |>
#         list_rbind()
#     t1 <- Sys.time()
#     print(t1 - t0)
#     return(out)
# }
#
# # Takes ~1 min
# # Note: do NOT attempt multithreading, just slows it down.
# # I tried both `future.apply` and `parallel::mclapply`.
# ea_info <- get_ea_info(ea_counts, ea_snapes)
#
#
# # This tells us that these alternative minor alleles are always the reference
# mean(ea_info$is_ref)
# # [1] 1
#
# # ... and that the allele frequency is always at or just under 0.5.
# min(ea_info$freq_min); max(ea_info$freq_max)
# # [1] 0.4385
# # [1] 0.5





