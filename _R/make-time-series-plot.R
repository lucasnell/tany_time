
#'
#' Make plot of candidate allele frequencies through time
#'

library(tidyverse)
library(ggtext)
library(patchwork)

theme_set(theme_classic() %+replace%
              theme(strip.background = element_blank(),
                    strip.text = element_markdown(size = 10),
                    axis.title = element_markdown(color = "black", size = 11),
                    axis.text = element_markdown(color = "black", size = 9),
                    axis.ticks = element_line(color = "black"),
                    legend.background = element_blank(),
                    plot.title = element_markdown(size = 14, hjust = 0.5)))


freq_df <- read_csv("_data/timeseries-data.csv", col_types = "fcificddddd") |>
    mutate(date = year + ifelse(gen == "a", 0, 0.5),
           logN = log1p(abund))

abund_df <- read_csv("_data/tany-abundances.csv", col_types = "icdd") |>
    mutate(date = year + ifelse(gen == "a", 0, 0.5),
           abund = N) |>
    select(year, gen, date, abund, logN) |>
    filter(date >= min(freq_df$date))

max_f <- max(freq_df$freq)
min_ln <- min(abund_df$logN)
max_ln <- max(abund_df$logN)

# logN --> freq scale
trans <- function(x) (x - min_ln) * max_f / (max_ln - min_ln)
# freq --> logN scale
itrans <- function(y) min_ln + y * (max_ln - min_ln) / max_f



frq_ts_p <- freq_df |>
    ggplot(aes(date, freq)) +
    geom_ribbon(data = abund_df |>
                    mutate(freq = trans(logN)),
                aes(ymin = 0, ymax = freq), fill = "gray80") +
    geom_line(aes(color = mrk)) +
    geom_point(aes(color = mrk)) +
    scale_y_continuous("Allele frequency", sec.axis = sec_axis(itrans, "log(1 + abundance)")) +
    scale_color_viridis_d(begin = 0.2, guide = "none") +
    theme(axis.title.y.right = element_text(color = "gray70"),
          axis.text.y.right = element_text(color = "gray70"),
          axis.ticks.y.right = element_line(color = "gray70"),
          axis.title.x = element_blank())
frq_ts_p
ggsave("~/Desktop/freq-abund-time.pdf", frq_ts_p, width = 6, height = 4)



freq_df |>
    ggplot(aes(logN, freq)) +
    geom_point(aes(color = mrk)) +
    stat_smooth(aes(color = mrk), formula = y ~ x, method = "lm", se = TRUE) +
    facet_wrap(~ mrk, nrow = 6, scales = "free_y") +
    theme(strip.text = element_blank()) +
    scale_color_viridis_d(begin = 0.2, guide = "none")


#' SNPs for two focal genes:
#'
#' 1. Q9VUD3 (sox21b), male genital, fertility
#' 2. Q9VDD2 (SNF4Agamma), responses to starvation
#'
focal_snps <- "~/Desktop/Tgraci-genomics/Tgraci-funct-Dmel.csv" |>
    read_csv(col_types = cols()) |>
    filter(grepl("Q9VUD3|Q9VDD2", saccver)) |>
    mutate(uniprot = str_split(saccver, "\\|") |> map_chr(\(x) x[[2]]),
           gene = qaccver) |>
    select(gene, uniprot) |>
    left_join(read_rds("_data/candi_funcs_df.rds") |>
                  group_by(gene) |>
                  filter(bf_db == max(bf_db)) |>
                  ungroup() |>
                  select(contig, pos, mrk, gene, dist),
              by = "gene") |>
    mutate(name = case_when(uniprot == "Q9VUD3" ~ "sox21b",
                            uniprot == "Q9VDD2" ~ "SNF4A&gamma;",
                            TRUE ~ NA))


focal_freq_p <- freq_df |>
    filter(mrk %in% focal_snps$mrk) |>
    mutate(name = map_chr(mrk, \(i) focal_snps$name[focal_snps$mrk == i])) |>
    ggplot(aes(logN, freq)) +
    geom_point(aes(color = mrk)) +
    stat_smooth(aes(color = mrk), formula = y ~ x, method = "lm", se = TRUE) +
    facet_wrap(~ name, ncol = 1) +
    theme(strip.text = element_markdown(size = 8, face = "italic"),
          axis.text = element_markdown(size = 6),
          axis.title = element_blank()) +
    scale_color_viridis_d(begin = 0.2, guide = "none")

# ggsave("_plots/focal-freqs.pdf", focal_freq_p, width = 3.65, height = 1.2)
# ggsave("_plots/focal-freqs.pdf", focal_freq_p, width = 1.5, height = 3)

