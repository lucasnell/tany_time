
#'
#' Summarize Kraken2 reports from `_bash/01-variant-calling/06-kraken`
#' (a random 1M unmapped reads per pool, classified against PlusPF-16).
#'
#' Kraken2 report columns: percent of reads in the clade, reads in the
#' clade, reads assigned directly to the taxon, rank code (U unclassified,
#' R root, D domain, K kingdom, P phylum, ..., G genus, S species, with
#' digits for intermediate ranks), NCBI taxon ID, and name (indented by
#' depth).
#'
#' Percentages are of the subsampled *unmapped* reads. If pool QC metrics
#' are available (`pool-qc-files.tsv` from `02-allele-freqs/03-pool-qc`),
#' the script also expresses them as a share of all reads:
#' percent / 100 * (1 - map_rate).
#'
#' Inputs: reports extracted from `kraken.tar` into `kraken_dir`.
#' Outputs: `_data/kraken-groups.csv` (key groups per pool),
#' `_data/kraken-top-genera.csv` (most common genera per pool), and
#' `_plots/kraken-groups.pdf`.
#'

library(tidyverse)

kraken_dir <- "~/_data/_tany-time/kraken"
qc_files <- "~/_data/_tany-time/pool-qc/pool-qc-files.tsv"

#' Groups to report, as (label, rank code, name). Fungi are a kingdom in
#' NCBI taxonomy; Homo sapiens is reported because handling contamination
#' is common.
groups <- tribble(
    ~label,          ~rank, ~name,
    "Unclassified",  "U",   "unclassified",
    "Bacteria",      "D",   "Bacteria",
    "Archaea",       "D",   "Archaea",
    "Viruses",       "D",   "Viruses",
    "Eukaryota",     "D",   "Eukaryota",
    "Fungi",         "K",   "Fungi",
    "Homo sapiens",  "S",   "Homo sapiens")

N_TOP <- 10  # genera per pool in the top-genera table



# =============================================================================*
# Read reports ----
# =============================================================================*

read_report <- function(fn) {
    read_tsv(fn, col_names = c("percent", "clade_reads", "direct_reads",
                               "rank", "taxid", "name"),
             col_types = "dddcic", trim_ws = FALSE, progress = FALSE) |>
        mutate(name = str_trim(name))
}

rep_df <- list.files(kraken_dir, "_kraken-report.txt$", full.names = TRUE) |>
    set_names(\(x) str_remove(basename(x), "_kraken-report.txt$")) |>
    map(read_report) |>
    list_rbind(names_to = "pool")

stopifnot(nrow(rep_df) > 0)

# Total reads classified per pool (unclassified + root):
totals <- rep_df |>
    filter(rank %in% c("U", "R")) |>
    group_by(pool) |>
    summarize(n_reads = sum(clade_reads), .groups = "drop")



# =============================================================================*
# Key groups ----
# =============================================================================*

groups_df <- crossing(pool = unique(rep_df$pool), groups) |>
    left_join(rep_df |> select(pool, rank, name, clade_reads),
              by = c("pool", "rank", "name")) |>
    mutate(clade_reads = replace_na(clade_reads, 0)) |>
    left_join(totals, by = "pool") |>
    mutate(percent = 100 * clade_reads / n_reads) |>
    select(pool, label, clade_reads, n_reads, percent)

if (file.exists(qc_files)) {
    map_df <- read_tsv(qc_files, col_types = "cddddddd", na = "NA") |>
        select(pool, map_rate)
    groups_df <- groups_df |>
        left_join(map_df, by = "pool") |>
        mutate(percent_of_all = percent * (1 - map_rate))
} else {
    message("No pool QC metrics found; reporting % of unmapped reads only")
}

write_csv(groups_df, "_data/kraken-groups.csv")

groups_df |>
    select(pool, label, percent) |>
    pivot_wider(names_from = label, values_from = percent) |>
    arrange(desc(Bacteria)) |>
    mutate(across(-pool, \(x) round(x, 2))) |>
    print(n = Inf, width = Inf)



# =============================================================================*
# Top genera ----
# =============================================================================*

top_df <- rep_df |>
    filter(rank == "G") |>
    left_join(totals, by = "pool") |>
    mutate(percent = 100 * clade_reads / n_reads) |>
    group_by(pool) |>
    slice_max(clade_reads, n = N_TOP, with_ties = FALSE) |>
    ungroup() |>
    select(pool, genus = name, taxid, clade_reads, percent)

write_csv(top_df, "_data/kraken-top-genera.csv")



# =============================================================================*
# Plot ----
# =============================================================================*

#' Domain-level composition of unmapped reads (Fungi and Homo sapiens are
#' within Eukaryota, so they're shown in a second panel).
p <- groups_df |>
    mutate(panel = ifelse(label %in% c("Fungi", "Homo sapiens"),
                          "Within Eukaryota", "Domains"),
           pool = fct_reorder(pool, percent * (label == "Bacteria"),
                              .fun = sum)) |>
    ggplot(aes(percent, pool, fill = label)) +
    geom_col() +
    facet_wrap(~ panel, scales = "free_x") +
    scale_fill_brewer(NULL, palette = "Dark2") +
    labs(x = "% of classified + unclassified unmapped reads", y = NULL) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 5))

ggsave("_plots/kraken-groups.pdf", p, width = 9, height = 8)
