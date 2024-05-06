source("renv/activate.R")

library(tidyverse)

raw_data <- read_csv(snakemake@input[[1]])
#raw_data <- read_csv("results/data/prepared/raw_abundance.csv")

# 1. Convert glycans-----
# Convert glycan strings from byonic format into condensed format.
# e.g. Hex(7)HexNAc(6)dHex(1)NeuAc[+13.0316](4) -> H7N6F1S4
# (The byonic format is from the result of GlyHunter.)
convert <- function(comp) {
  comp = str_replace(comp, "\\[.*?\\]", "")
  comp = str_replace(comp, "Hex\\((\\d+)\\)", "H\\1")
  comp = str_replace(comp, "HexNAc\\((\\d+)\\)", "N\\1")
  comp = str_replace(comp, "dHex\\((\\d+)\\)", "F\\1")
  comp = str_replace(comp, "NeuAc\\((\\d+)\\)", "S\\1")
  return(comp)
}

converted <- raw_data |> 
  mutate(glycan = convert(glycan))

# 2. Filter glycans-----
# Remove glycans with missing values in more than 50% of samples.
filtered_1 <- converted |> 
  group_by(glycan) |> 
  mutate(missing_prop = mean(is.na(value))) |> 
  filter(missing_prop < 0.5) |> 
  ungroup() |> 
  select(-missing_prop)

n_distinct(converted$glycan)
n_distinct(filtered_1$glycan)
# Before filtering: 72 glycans
# After filtering: 62 glycans

# 3. Filter samples-----
to_delete_mannual <- c("D231", "D219", "D243", "D194", "D212")  # Outlier samples from manual inspection
to_delete_missing <- filtered_1 |>   # Outlier samples with too many missing values
  summarise(missing_prop = mean(is.na(value)), .by = sample) |> 
  filter(missing_prop >= 0.5) |> 
  pull(sample)
to_delete <- c(to_delete_mannual, to_delete_missing)
filtered_2 <- filtered_1 |> 
  filter(!sample %in% to_delete)

# 4. Normalize abundance-----
# Normalize abundance values using median quotient normalization.
normalized <- filtered_2 |> 
  # Replace missing values with 0
  # This step is essential for a robust normalization based on median quotient.
  # If missing values are not replaced, the medians will be based on different sets of samples.
  # This will result in a biased normalization.
  mutate(
    missing = is.na(value),
    value = replace_na(value, 0)
  ) |> 
  rename(gid = sample) |>  # This is required by glycanr
  glycanr::medianquotientnorm() |> 
  rename(sample = gid) |>   # Rename back to sample
  mutate(value = if_else(missing, NA, value)) |> 
  select(-missing)

# 5. Impute glycans-----
wide_normalized <- normalized |> 
  pivot_wider(names_from = glycan, values_from = value)
imputed <- wide_normalized |> 
  select(-sample) |> 
  # Impute with KNN
  VIM::kNN(k = 10) |> as_tibble() |> 
  select(-ends_with("imp")) |> 
  mutate(sample = wide_normalized$sample) |> 
  pivot_longer(-sample, names_to = "glycan", values_to = "value")

write_csv(imputed, snakemake@output[[1]])