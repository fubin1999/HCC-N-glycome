source("renv/activate.R")

library(tidyverse)

raw_data <- read_csv(snakemake@input[[1]])
# raw_data <- read_csv("results/data/prepared/raw_abundance.csv")

# 1. Convert to long-----
long_data <- raw_data |>
  pivot_longer(-sample, names_to = "glycan", values_to = "value")

# 2. Filter samples-----
# Filter outlier samples based on the number of glycans detected.
to_delete_1 <- c("S231", "S219", "S243", "S194", "S212")
to_delete_2 <- long_data |>
  summarise(na_prop = mean(is.na(value)), .by = sample) |> 
  filter(na_prop >= 0.5) |> 
  pull(sample)
to_delete <- c(to_delete_1, to_delete_2)

filtered_1 <- long_data |>
  filter(!sample %in% to_delete)

# 3. Filter glycans-----
# Remove glycans with missing values in more than 50% of samples.
filtered_2 <- filtered_1 |> 
  group_by(glycan) |> 
  mutate(missing_prop = mean(is.na(value))) |> 
  filter(missing_prop < 0.5) |> 
  ungroup() |> 
  select(-missing_prop)

# Before filtering: 72 glycans
# After filtering: 62 glycans

# 4. Impute glycans-----
# Impute missing values using the the minimum value of each sample / 2.
imputed <- filtered_2 |> 
  group_by(sample) |> 
  mutate(min_value = min(value, na.rm = TRUE)) |> 
  ungroup() |> 
  mutate(value = if_else(is.na(value), min_value / 2, value)) |> 
  select(-min_value)

# 5. Normalize abundance-----
# Normalize abundance values using median quotient normalization.
normalized <- imputed |> 
  rename(gid = sample) |>  # This is required by glycanr::
  glycanr::medianquotientnorm() |> 
  rename(sample = gid)  # Rename back to sample

# 6. Convert back to wide-----
final_data <- normalized |>
  pivot_wider(names_from = glycan, values_from = value)

# Save result-----
write_csv(final_data, snakemake@output[[1]])