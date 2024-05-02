library(tidyverse)

files <- unlist(snakemake@input)
plates <- read_csv(files)
groups <- plates |> 
    mutate(group = if_else(
        raw == "QC", "QC", str_sub(raw, 1, 1)
    )) |> 
    select(sample, group)
write_csv(groups, snakemake@output[[1]])