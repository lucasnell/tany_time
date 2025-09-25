
#'
#' Make dataset of allele frequencies for snps associated with abundance,
#' for use in time series analysis aimed at looking for balancing selection.
#'

library(tidyverse)
library(Rcpp)


sourceCpp("_R/snape-helpers.cpp")


snape_file <- function(x){
    paste0("~/_data/_tany-time/snape/snape_masked_time/", x)
}
count_file <- function(x){
    paste0("~/_data/_tany-time/snape/count_masked_time/", x)
}

samp_df <- read_csv("_data/full-sample-info.csv", col_types = "cDccidd")

pools <- read_lines(snape_file("snape_masked_time.names.gz"))
n_adults <- pools |>
    set_names() |>
    map_int(\(x) samp_df$n_adults[samp_df$biotech_id == x]) |>
    as.list()

snape_df <- read_tsv(snape_file("snape_masked_time_noblanks.sync.gz"),
                     col_names = c("contig", "pos", "ref", pools),
                     col_types = paste(c("cic", rep("c", length(pools))),
                                       collapse = ""),
                     progress = FALSE) |>
    mutate(contig = factor(contig, levels = paste0("CONTIG", 1:45)),
           across(all_of(pools), get_snape_af))

count_df <- read_tsv(count_file("count_masked_time_noblanks.sync.gz"),
                     col_names = c("contig", "pos", "ref", pools),
                     col_types = paste(c("cic", rep("c", length(pools))),
                                       collapse = ""),
                     progress = FALSE) |>
    mutate(contig = factor(contig, levels = paste0("CONTIG", 1:45)),
           across(all_of(pools), split_count_strings))

freq_df <- snape_df |>
    pivot_longer(all_of(pools), names_to = "pool", values_to = "freq") |>
    full_join(count_df |>
                  mutate(across(all_of(pools), \(x) map_int(x, sum))) |>
                  pivot_longer(all_of(pools), names_to = "pool",
                               values_to = "n_reads"),
              by = c("contig", "pos", "ref", "pool")) |>
    mutate(n_haps = 2L * map_int(pool, \(p) n_adults[[p]]),
           n_eff = (n_reads * n_haps - 1) / (n_reads + n_haps))


# Candidate SNPs from BayPass:
candi_df <- read_rds("_data/candi_df.rds")

candi_locs <- paste(candi_df$contig, candi_df$pos, sep = "_")

candi_freqs <- freq_df |>
    mutate(location = paste(contig, pos, sep = "_")) |>
    filter(location %in% candi_locs) |>
    select(-location)

# Should return TRUE
isTRUE(length(pools) == (candi_freqs |>
                             group_by(contig, pos) |>
                             summarize(n = n(), .groups = "drop") |>
                             getElement("n") |>
                             unique()))


tany_pop_df <- read_csv("_data/tany-abundances.csv", col_types = "icdd")


ts_df <- candi_freqs |>
    mutate(mrk = map2_int(contig, pos, \(c, p) {
        candi_df$mrk[candi_df$contig == c & candi_df$pos == p]
    }),
           info = pool |> str_remove_all("_S.*") |> str_split("-"),
           site = map_chr(info, \(x) x[[1]]),
           year = map_int(info, \(x) {
               y <- as.integer(x[[3]])
               if (y > 70) y <- y + 1900L
               else y <- y + 2000L
               return(y)
           }),
           gen = map_chr(info, \(x) letters[as.integer(x[[2]])]),
    abund = map2_dbl(year, gen, \(y, g) {
        tany_pop_df$N[tany_pop_df$year == y & tany_pop_df$gen == g]
    })) |>
    select(mrk, contig, pos, pool, year, gen, freq:n_eff, abund)


#' Column descriptions:
#'     1. mrk - Number for this locus, used in baypass
#'     2. contig - Locus contig name
#'     3. pos - Locus positions on contig
#'     4. pool - Name of sampling pool
#'     5. year - Year sample was collected
#'     6. gen - Generation sample was collected (either a or b)
#'     7. freq - Minor allele frequency at this locus for this sample
#'     8. n_reads - Total number of mapped reads for this sample at this locus
#'     9. n_haps - Number of haploid genomes (2 * number of adults) in the sample
#'     10. n_eff - Effective sample size ((n_reads * n_haps - 1) / (n_reads + n_haps))
#'     11. abund - Tanytarsus gracilentus abundance in the generation coinciding with this sample
#'

write_csv(ts_df, "_data/timeseries-data.csv")


