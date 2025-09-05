
library(tidyverse)


#'
#' Get sample names copied to CSV file from HTML file provided by biotech center
#'
#' Simple names are those that show up in data I sent to biotech center
#' (in CSV files below).
#'
samp_names <- read_csv("_data/biotech_ids.csv", col_types = "c",
                       col_names = "name") |>
    mutate(simple = str_remove(name, "_S.*")) |>
    (\(x) {
        z <- as.list(x$name)
        names(z) <- x$simple
        return(z)
    })()

#'
#' Emergence trap locations:
#'
trap_df <- read_csv("_data/trap_locations.csv", col_types = cols()) |>
    # Simplify site names:
    mutate(trap = case_when(trap == "Syðri Neslönd" ~ "SN",
                            trap == "Kálfaströnd" ~ "KS",
                            trap == "Haganes" ~ "H",
                            trap == "Vindbelgur" ~ "V",
                            TRUE ~ NA_character_))


#' Archived tany info:
archive_df <- read_csv("_data/tany-from-archive.csv", col_types = cols()) |>
    mutate(date = as.Date(sprintf("%i-%i-%i", col_year, col_month, col_day)),
           lake = "Mývatn") |>
    filter(!is.na(biotech_id)) |>
    filter(biotech_id %in% names(samp_names)) |>
    #'
    #' These samples were prepped together:
    #'
    #' KS-2-90 and KS-2-90-B
    #' SN-2-91 and SN-2-91-B
    #' KS-2-83 and KS-2-83-B
    #'
    mutate(biotech_id = ifelse(biotech_id %in% c("KS-2-90-B", "SN-2-91-B", "KS-2-83-B"),
                               str_remove(biotech_id, "-B"),
                               biotech_id)) |>
    group_by(biotech_id, date, lake, site) |>
    summarize(n_adults = as.integer(sum(num_males)), .groups = "drop") |>
    mutate(lat = map_dbl(site, ~ trap_df$lat[trap_df$trap == .x]),
           lon = map_dbl(site, ~ trap_df$lon[trap_df$trap == .x]))

# Should be FALSE
any(is.na(archive_df))



spat_df <- read_csv("_data/tany-other-lakes.csv", col_types = cols()) |>
    mutate(date = as.Date(sprintf("2019-%i-%i", col_month, col_day)),
           across(col_month:n_tany, as.integer)) |>
    rename(month = col_month, day = col_day, n_adults = n_tany) |>
    filter(!is.na(biotech_id)) |>
    filter(biotech_id %in% names(samp_names)) |>
    select(biotech_id, date, lake, sample, month, day, n_adults) |>
    #'
    #' Had to add this bc in `tany-other-lakes-gps.csv`, before June 20th,
    #' when only a single sample was taken from a site on a particular day,
    #' I called it sample # 1.
    #' After June 20th, I just left it blank.
    #' In `tany-other-lakes.csv`, I just always left these blank.
    #'
    mutate(sample = ifelse(is.na(sample) & month == 6 & day < 20, 1, sample)) |>
    left_join(read_csv("_data/tany-other-lakes-gps.csv", col_types = cols()) |>
                  select(lake, sample, month, day, lat, lon) |>
                  mutate(month = map_int(month, \(x) which(month.name == x)),
                         day = as.integer(day)),
              by = c("lake", "sample", "month", "day")) |>
    mutate(site = map_chr(str_split(lake, " - "), ~ tail(.x, 1)),
           lake = map_chr(str_split(lake, " - "), ~ .x[[1]]),
           site = case_when(lake != "Mývatn" ~ NA_character_,
                            site == "Mývatn" ~ "BR",
                            TRUE ~ site)) |>
    # "Blonduos" is the lake's nearest town; the lake is "Grafarvatn"
    mutate(lake = ifelse(lake == "Blonduos", "Grafarvatn", lake)) |>
    select(biotech_id, date, lake, site, n_adults, lat, lon)



#'
#' Combine and make biotech IDs match with downstream naming:
#'
full_df <- bind_rows(archive_df, spat_df) |>
    mutate(biotech_id = map_chr(biotech_id, \(x) samp_names[[x]]))



write_csv(full_df, "_data/full-sample-info.csv")

