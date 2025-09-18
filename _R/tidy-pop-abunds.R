

library(tidyverse)


tany_pop_df <- read_csv("_data/chironomid-abundances.csv", col_types = "iddddd") |>
    select(year, starts_with("tany")) |>
    pivot_longer(starts_with("tany"), names_to = "gen", values_to = "N") |>
    mutate(gen = str_remove_all(gen, "^tany_") |> factor(),
           logN = log1p(N))

write_csv(tany_pop_df, "_data/tany-abundances.csv")

