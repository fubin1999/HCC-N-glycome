source("renv/activate.R")

library(tidyverse)

plates <- read_csv(snakemake@input[[1]])
groups <- plates |> 
    mutate(group = if_else(
        raw_sample == "QC", "QC", str_sub(raw_sample, 1, 1)
    )) |> 
    select(sample, group)
write_csv(groups, snakemake@output[[1]])