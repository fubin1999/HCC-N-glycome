library(tidyverse)

# processed_file <- "results/data/prepared/processed_abundance.csv"
processed_file <- snakemake@input[[1]]
data <- read_csv(processed_file) |> 
  pivot_wider(names_from = "glycan", values_from = "value")

glycans <- setdiff(colnames(data), "sample")

result <- data |> 
  mutate(Gtest = log((H6N5F1S1 + H6N5F1S2 + H6N5F1S3) / (H5N5F1 + H6N5F1S1 + H6N5F1S2 + H6N5F1S3))) |> 
  select(sample, Gtest)

write_csv(result, snakemake@output[[1]])