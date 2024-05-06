# Combine processed abundance data with groups, and filter out QC samples.

source("renv/activate.R")

library(tidyverse)

abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

result <- groups |>
  left_join(abundance, by = "sample") |> 
  filter(group != "QC") |> 
  mutate(group = factor(group, levels = c("H", "M", "Y", "C")))

write_rds(result, snakemake@output[[1]])