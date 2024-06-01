library(tidyverse)

pos_file <- snakemake@input[[1]]
glyhunter_result_file <- file.path(snakemake@input[[2]], "summary_area.csv")

plate_no <- as.numeric(str_extract(glyhunter_result_file, "plate(\\d+)", group = 1))

pos_df <- read_csv(pos_file) |> 
  filter(plate == plate_no) |> 
  select(-plate)
glyhunter_df <- read_csv(glyhunter_result_file)

glyhunter_long <- glyhunter_df |>
  pivot_longer(-glycan, names_to="sample", values_to="value")

result <- glyhunter_long |>
  separate_wider_regex(sample, c(".*_", position = ".*", "_1")) |>
  left_join(pos_df, by = "position") |>
  select(sample, glycan, value)

write_csv(result, snakemake@output[[1]])