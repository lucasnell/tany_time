
#'
#' Summarize Kraken2 reports from `_bash/01-variant-calling/05-kraken`:
#' - PlusPF-16 (`<pool>_kraken-report.txt`): a random 1M unmapped reads per
#'   pool; overall composition (bacteria, fungi, human, viruses).
#' - SILVA SSU rRNA (`<pool>_silva-report.txt`): all unmapped reads;
#'   eukaryotic parasites (last section).
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
#' `_data/kraken-top-genera.csv` (most common genera per pool),
#' `_plots/kraken-groups.pdf`, `_data/kraken-silva-groups.csv` (parasite
#' groups per pool), and `_data/kraken-silva-top-eukaryotes.csv`.
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

#' Names are indented by 2 spaces per level; `depth` keeps that level so
#' lineages can be rebuilt (see `in_lineage`).
read_report <- function(fn) {
    read_tsv(fn, col_names = c("percent", "clade_reads", "direct_reads",
                               "rank", "taxid", "name"),
             col_types = "dddcic", trim_ws = FALSE, progress = FALSE) |>
        mutate(depth = (nchar(name) - nchar(str_trim(name, "left"))) %/% 2,
               name = str_trim(name))
}

#' For rows of one report (in file order), is `target` the taxon itself or
#' one of its ancestors?
in_lineage <- function(name, depth, target) {
    stack <- character(0)
    out <- logical(length(name))
    for (i in seq_along(name)) {
        stack <- c(stack[seq_len(depth[i])], name[i])
        out[i] <- target %in% stack
    }
    out
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



# =============================================================================*
# SILVA: eukaryotic parasites ----
# =============================================================================*

#' SILVA taxonomy names differ from NCBI's, so each group lists candidate
#' names; reads are summed over the names found within a group (clade
#' counts, so a group's names must not be nested in each other). Groups
#' can overlap (e.g., "Fungi (all)" includes Microsporidia and
#' Blastocladiomycota in SILVA), so don't sum across groups.
#' The script reports which names matched nothing in any report, so
#' unmatched groups aren't mistaken for absences.
#'
#' Host rRNA (Chironomidae) is the denominator: parasite reads per host
#' rRNA read compare infection signal among pools with different read
#' numbers. Parasite DNA may come from other specimens in the same trap
#' preservative, so detection is not proof that T. gracilentus was infected.
silva_groups <- tribble(
    ~label,                ~names,
    "Chironomidae (host)", list("Chironomidae"),
    "Nematoda",            list("Nematoda"),
    "Mites",               list("Acari", "Acariformes", "Parasitiformes"),
    "Microsporidia",       list("Microsporidia"),
    "Blastocladiomycota",  list("Blastocladiomycota", "Blastocladiomycetes"),
    "Entomophthorales",    list("Entomophthorales", "Entomophthoromycotina"),
    "Oomycetes",           list("Peronosporomycetes", "Oomycota"),
    "Fungi (all)",         list("Fungi"),
    "Platyhelminthes",     list("Platyhelminthes"),
    "Apicomplexa",         list("Apicomplexa")) |>
    unnest_longer(names) |>
    mutate(names = unlist(names))

silva_df <- list.files(kraken_dir, "_silva-report.txt$", full.names = TRUE) |>
    set_names(\(x) str_remove(basename(x), "_silva-report.txt$")) |>
    map(read_report) |>
    list_rbind(names_to = "pool")

if (nrow(silva_df) == 0) {
    message("No SILVA reports found")
} else {

    unmatched <- setdiff(silva_groups$names, silva_df$name)
    if (length(unmatched) > 0) {
        message("SILVA names not found in any report: ",
                paste(unmatched, collapse = ", "))
    }

    silva_tot <- silva_df |>
        filter(rank %in% c("U", "R")) |>
        group_by(pool) |>
        summarize(n_reads = sum(clade_reads), .groups = "drop")

    silva_groups_df <- crossing(pool = unique(silva_df$pool), silva_groups) |>
        left_join(silva_df |> select(pool, name, clade_reads) |>
                      distinct(pool, name, .keep_all = TRUE),
                  by = c("pool", "names" = "name")) |>
        mutate(clade_reads = replace_na(clade_reads, 0)) |>
        group_by(pool, label) |>
        summarize(clade_reads = sum(clade_reads), .groups = "drop") |>
        left_join(silva_tot, by = "pool")

    host_df <- silva_groups_df |>
        filter(label == "Chironomidae (host)") |>
        select(pool, host_reads = clade_reads)

    silva_groups_df <- silva_groups_df |>
        left_join(host_df, by = "pool") |>
        mutate(per_host_read = ifelse(host_reads > 0,
                                      clade_reads / host_reads, NA_real_))

    write_csv(silva_groups_df, "_data/kraken-silva-groups.csv")

    silva_groups_df |>
        filter(label != "Chironomidae (host)") |>
        select(pool, label, clade_reads) |>
        pivot_wider(names_from = label, values_from = clade_reads) |>
        print(n = Inf, width = Inf)

    #' Most common eukaryotic genera outside Insecta, per pool (lineages
    #' rebuilt from the report's indentation).
    silva_top <- silva_df |>
        group_by(pool) |>
        mutate(euk = in_lineage(name, depth, "Eukaryota"),
               insect = in_lineage(name, depth, "Insecta")) |>
        ungroup() |>
        filter(rank == "G", euk, !insect) |>
        left_join(silva_tot, by = "pool") |>
        left_join(host_df, by = "pool") |>
        mutate(per_host_read = ifelse(host_reads > 0,
                                      clade_reads / host_reads, NA_real_)) |>
        group_by(pool) |>
        slice_max(clade_reads, n = N_TOP, with_ties = FALSE) |>
        ungroup() |>
        select(pool, genus = name, taxid, clade_reads, per_host_read)

    write_csv(silva_top, "_data/kraken-silva-top-eukaryotes.csv")
}
