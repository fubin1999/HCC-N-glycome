source("renv/activate.R")

library(tidyverse)

files <- list.files(snakemake@input[[1]], full.names = TRUE, pattern = ".csv")
df <- read_csv(files, id = "plate") |> 
  rename(raw_sample = raw) |> 
  separate_wider_regex(plate, c(".*plate", plate = "\\d+", "\\.csv"))
write_csv(df, snakemake@output[[1]])