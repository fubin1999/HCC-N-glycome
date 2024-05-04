source("renv/activate.R")

library(tidyverse)

plates <- read_csv(snakemake@input[["plates"]])
abundance <- read_csv(snakemake@input[["abundance"]])
groups <- plates |> 
    mutate(group = if_else(
        raw_sample == "QC", "QC", str_sub(raw_sample, 1, 1)
    )) |> 
    select(sample, group) |> 
    semi_join(abundance, by = "sample")
write_csv(groups, snakemake@output[[1]])