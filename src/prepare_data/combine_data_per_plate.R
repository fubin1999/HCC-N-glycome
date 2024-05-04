source("renv/activate.R")
library(tidyverse)

files <- unlist(snakemake@input)
combined <- read_csv(files) |> 
  complete(sample, glycan, fill = list(value = NA))
write_csv(combined, snakemake@output[[1]])