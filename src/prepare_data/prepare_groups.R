library(tidyverse)

plates <- read_csv(snakemake@input[["plates"]])

groups <- plates |>
    mutate(group = if_else(
        raw_sample == "QC", "QC", str_sub(raw_sample, 1, 1)
    )) |> 
    mutate(group = case_match(group, "H" ~ "HC", "M" ~ "CHB", "Y" ~ "LC", "C" ~ "HCC", "QC" ~ "QC")) |> 
    select(sample, group)

write_csv(groups, snakemake@output[[1]])