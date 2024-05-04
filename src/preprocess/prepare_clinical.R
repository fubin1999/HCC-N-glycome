source("renv/activate.R")

library(tidyverse)

raw_clinical <- read_csv(snakemake@input[["clinical"]])
abundance <- read_csv(snakemake@input[["abundance"]])
plates <- read_csv(snakemake@input[["plates"]])

# Convert sample names and keep only samples that have abundance data
clinical <- raw_clinical |> 
  left_join(plates |> select(raw_sample, sample), by = "raw_sample") |> 
  select(sample, everything()) |>
  select(-raw_sample) |> 
  semi_join(abundance |> distinct(sample), by = "sample") |> 
  separate_wider_regex(sample, c("D", sample_no = "\\d+"), cols_remove = FALSE) |> 
  mutate(sample_no = as.integer(sample_no)) |> 
  arrange(sample_no) |> 
  select(-sample_no)

# Clean up clinical data
clean_clinical <- clinical |> 
  mutate(across(-c(sample, sex), ~ parse_number(as.character(.x)))) |> 
  mutate(across(-c(sample, sex), ~ replace_na(.x, 0))) |> 
  distinct(sample, .keep_all = TRUE)

write_csv(clean_clinical, snakemake@output[[1]])
