
#'
#' ```bash
#' # Beforehand:
#' # Download tar files from `01-fastp` step into the `~/_data/trimmed_reports`
#' # directory.
#' # Install `multiqc` into new conda environment (v 1.30):
#'
#' conda create -n multiqc-env multiqc=1.30
#'
#' conda activate multiqc-env
#'
#'
#' # Go to this directory of reports and un-tar all files
#' cd ~/_data/trimmed_reports
#' for f in *.tar.gz; do tar -xzf $f; done
#' rm *.tar.gz
#'
#'
#' # Adjust names for fastp reports so multiqc can parse sample names:
#' for f in trimmed_*; do
#'     sample=${f/trimmed_/}
#'     mv $f $sample
#'     cd $sample
#'     mv fastp.html ${sample}.html
#'     mv fastp.json ${sample}.json
#'     cd ..
#' done
#'
#' # Run MultiQC
#' multiqc .
#'
#' # If no issues above, open the report
#' open multiqc_report.html
#' ```
#'


library(tidyverse)
library(readxl)
library(ggtext)
library(patchwork)

theme_set(theme_classic() %+replace%
              theme(strip.background = element_blank(),
                    strip.text = element_markdown(size = 10),
                    axis.title = element_markdown(color = "black", size = 11),
                    axis.text = element_markdown(color = "black", size = 9),
                    axis.ticks = element_line(color = "black"),
                    legend.background = element_blank(),
                    plot.title = element_markdown(size = 14, hjust = 0.5),
                    panel.grid.major = element_line()))



# read_tsv("~/_data/trimmed_reports/multiqc_data/multiqc_fastqc.txt",
#          col_types = cols()) |>
#     filter()


if (!file.exists("_data/bed-summs.rds")) {
    if (!dir.exists("/Volumes/Lucas512")) stop("Insert flash drive!")
    # Takes ~ 30 min
    bed_summs <- "/Volumes/Lucas512/other/reads/illumina/dna/bwa" |>
        list.files(".bed.gz$", full.names = TRUE) |>
        map(\(x) {
            x_bed_summ <- x |>
                read_tsv(col_names = c("contig", "start", "end", "coverage"),
                         col_types = "ciii", progress = FALSE) |>
                mutate(n_bases = end - start) |>
                summarize(p0 = sum(n_bases[coverage == 0]) / sum(n_bases),
                          min = min(coverage),
                          max = max(coverage),
                          q20 = quantile(coverage, 0.2),
                          q80 = quantile(coverage, 0.8),
                          median = median(rep(coverage, n_bases)),
                          mean = sum(coverage * n_bases) / sum(n_bases),
                          log_mean = sum(log1p(coverage) * n_bases) / sum(n_bases)) |>
                mutate(sample = basename(x) |> str_remove("_bwa.*")) |>
                select(sample, everything())
            return(x_bed_summ)
        },
        .progress = list(clear = FALSE,
                         format = paste("{cli::pb_bar} {cli::pb_percent}",
                                        "[{cli::pb_elapsed}] | ETA: {cli::pb_eta}"))) |>
        list_rbind()
    write_rds(bed_summs, "_data/bed-summs.rds")
} else {
    bed_summs <- read_rds("_data/bed-summs.rds")
}






# fastqc <- read_excel("~/Desktop/multiqc_fastqc.xlsx") |>
#     select(Sample, `%GC`, total_deduplicated_percentage,
#            per_sequence_gc_content, sequence_duplication_levels,
#            overrepresented_sequences, adapter_content,
#            Total_bases_Mb, Coverage, to_use) |>
#     rename_with(\(x) tolower(x) |> str_remove("%")) |>
#     rename(flag_gc = per_sequence_gc_content,
#            flag_dup = sequence_duplication_levels,
#            flag_over = overrepresented_sequences,
#            flag_adapt = adapter_content) |>
#     mutate(dup = 100 - total_deduplicated_percentage,
#            read = str_sub(sample, -5, -5) |> as.integer(),
#            sample = str_remove_all(sample, "trimmed_|_L002.*")) |>
#     # select(read, sample, everything()) |>
#     # mutate(test = ifelse(flag_gc == "fail" | flag_dup == "fail", 0, 1))
#     group_by(sample) |>
#     summarize(gc = mean(gc), dup = mean(dup), coverage = sum(coverage),
#               across(starts_with("flag_"), \(x) {
#                   sum(x == "fail")
#               }),
#               to_use = mean(to_use)) |>
#     rename_with(\(x) gsub("flag_", "fails_", x)) |>
#     mutate(test = as.integer(fails_gc == 0 & fails_dup == 0 & gc < 40),
#            test2 = as.integer(gc < 50))

fastqc <- read_tsv("~/_data/trimmed_reports/multiqc_data/multiqc_fastqc.txt",
                   col_types = cols()) |>
    select(Sample, `%GC`, total_deduplicated_percentage,
           per_sequence_gc_content, sequence_duplication_levels,
           overrepresented_sequences, adapter_content) |>
    rename_with(\(x) tolower(x) |> str_remove("%")) |>
    rename(flag_gc = per_sequence_gc_content,
           flag_dup = sequence_duplication_levels,
           flag_over = overrepresented_sequences,
           flag_adapt = adapter_content) |>
    mutate(dup = 100 - total_deduplicated_percentage,
           read = str_sub(sample, -5, -5) |> as.integer(),
           sample = str_remove_all(sample, "trimmed_|_L002.*")) |>
    group_by(sample) |>
    summarize(gc = mean(gc), dup = mean(dup),
              across(starts_with("flag_"), \(x) {
                  sum(x == "fail")
              })) |>
    rename_with(\(x) gsub("flag_", "fails_", x))




qc_df <- full_join(fastqc, bed_summs, by = "sample")
qc_df



qc_df |>
    ggplot(aes(gc, log10(p0), color = factor(fails_gc))) +
    # Known GC content of assembly:
    geom_vline(xintercept = 32.5, color = "gray70", linewidth = 0.75) +
    geom_point() +
    xlab("GC percent") +
    scale_y_continuous("Proportion genome missing",
                       breaks = log10(c(0.025, 0.05, 0.1, 0.2, 0.4, 0.8)),
                       labels = c(0.025, 0.05, 0.1, 0.2, 0.4, 0.8)) +
    scale_color_brewer(palette = "Dark2")


qc_df |>
    ggplot(aes(gc, median, color = factor(fails_over))) +
    # Known GC content of assembly:
    geom_vline(xintercept = 32.5, color = "gray70", linewidth = 0.75) +
    geom_point() +
    xlab("GC percent") +
    scale_y_continuous("Median coverage") +
    scale_color_brewer(palette = "Dark2")


qc_df |>
    mutate(gc_over = factor(gc > 40)) |>
    ggplot(aes(log10(p0), median, color = factor(fails_gc))) +
    geom_point(aes(shape = gc_over)) +
    scale_x_continuous("Proportion genome missing",
                       breaks = log10(c(0.025, 0.05, 0.1, 0.2, 0.4, 0.8)),
                       labels = c(0.025, 0.05, 0.1, 0.2, 0.4, 0.8)) +
    scale_y_continuous("Median coverage") +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual("GC > 40%", values = c(19, 1))


# -------------------------------*
# Final samples to use ----
# -------------------------------*

# For time series:
qc_df |> filter(median >= 30, str_starts(sample, "SN-")) |>
    getElement("sample") |>
    paste(collapse = "\n") |> cat()

# SN-1-02_S55
# SN-1-07_S23
# SN-1-08_S24
# SN-1-09_S25
# SN-1-10_S26
# SN-1-13_S56
# SN-1-14_S27
# SN-1-97_S61
# SN-1-99_S28
# SN-2-00_S29
# SN-2-02_S30
# SN-2-07_S31
# SN-2-08_S32
# SN-2-15_S34
# SN-2-92_S35
# SN-2-97_S37
# SN-2-99_S38


# For space:
qc_df |> filter(median >= 30, !str_starts(sample, "KS-|H-|SN-")) |>
    getElement("sample") |>
    paste(collapse = "\n") |> cat()

# Ash-19_S5
# Blik-19_S6
# Blo-19_S7
# Ellid-19_S8
# Hrisatjorn_S11
# Lei-19_S13
# Lys-19_S14
# Mik-19_S15
# MikW-19_S16
# MyBR-19_S17
# MyKS-19-B_S18
# MySN-19_S19
# Ola-19_S20
# Pernuvatn_S21
# Skj-19_S22
# Sva-19_S39
# Tjornin_S40
# Vik-19_S41





# =============================================================================*
# =============================================================================*
# OLD CODE ----
# =============================================================================*
# =============================================================================*


# # Heatmap of pass / fail / warning output from multiqc report:
# heatmap <- read_tsv("~/_data/trimmed_reports/multiqc_data/fastqc-status-check-heatmap.txt",
#                     col_types = cols()) |>
#     rename(sample = `.`) |>
#     pivot_longer(-sample, names_to = "stat", values_to = "flag") |>
#     mutate(flag = case_when(flag == 1 ~ "pass",
#                             flag == 0.5 ~ "warn",
#                             flag == 0.25 ~ "fail",
#                             TRUE ~ NA_character_),
#            read = str_sub(sample, -5, -5) |> as.integer(),
#            sample = str_remove_all(sample, "trimmed_|_L002.*")) |>
#     mutate(stat = case_when(stat == "Sequence Duplication Levels" ~ "duplicate_flag",
#                             stat == "Overrepresented Sequences" ~ "overrep_flag",
#                             stat == "Adapter Content" ~ "adapter_flag",
#                             TRUE ~ stat)) |>
#     filter(stat %in% c("duplicate_flag", "overrep_flag", "adapter_flag")) |>
#     pivot_wider(names_from = stat, values_from = flag)
#
#
# heatmap |>
#     filter()




# Main table of fastq info:
stat_df <- read_excel("~/Box Sync/midges/fastq-stats.xlsx", sheet = "trimmed") |>
    rename(sample = `Sample Name`)

# They *should* be the same order, but better to check:
identical(stat_df$sample, heatmap$sample)


# write_csv(left_join(stat_df, heatmap, by = "sample"),
#           "~/Box Sync/midges/fastq-stats-with-flags.csv")



#' ```bash
#' export OUT_CSV=mpileup_summs.csv
#' for f in *_mpileup.txt.gz; do
#'     echo -n ${f} "started..."
#'     echo -n ${f/_mpileup.txt.gz/}"," >> ${OUT_CSV}
#'     gunzip -c ${f} | cut -f 4 | awk '{s+=$1} END {printf "%.0f\n", s}' \
#'         >> ${OUT_CSV}
#'     echo " and finished"
#' done
#' ```


# Total sequencing output across entire genome based on mpileup files
final_out_df <- read_csv("~/_data/final_coverage/mpileup_summs.csv",
                         col_names = c("biotech_id", "total_seq"),
                         col_types = "cd")

# Add biotech_id col to match
stat_df <- stat_df |>
    mutate(biotech_id = str_remove_all(sample, "trimmed_|_L002_R2_001|_L002_R1_001"))

identical(sort(unique(stat_df$biotech_id)), final_out_df$biotech_id)

# Dividing by two bc the `stat_df` is split by FASTQ file (2 per sample),
# while `final_out_df` is split by BAM file (1 per sample)
stat_final_out_df <- left_join(stat_df,
                               mutate(final_out_df, total_seq = total_seq / 2),
                               by = "biotech_id")
# write_csv(stat_final_out_df,
#           "~/Box Sync/midges/fastq-stats-with-final-out.csv")


# Now combine these in the excel file.
# To see which are to be analyzed:

stat_df <- read_excel("~/Box Sync/midges/fastq-stats.xlsx", sheet = "trimmed") |>
    rename(sample = `Sample Name`) |>
    mutate(biotech_id = str_remove_all(sample, "trimmed_|_L002_R2_001|_L002_R1_001")) |>
    filter(to_use == 1)

# Make sure there aren't any loners:
stat_df |>
    group_by(biotech_id) |>
    summarize(N = n(), .groups = "drop") |>
    .[["N"]] |> unique()

an_names <- stat_df |>
    .[["biotech_id"]] |>
    unique()

temp_filt <- Reduce(`|`, lapply(c("^SN-", "^KS-", "^H-"), grepl, x = an_names))
temp_samples <- an_names[temp_filt]
paste0(temp_samples, "*", collapse = " ")


