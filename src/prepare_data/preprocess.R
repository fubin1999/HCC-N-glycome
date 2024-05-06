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

# 2. Filter samples-----
# Filter outlier samples based on the number of glycans detected.
to_delete <- c("D231", "D219", "D243", "D194", "D212")
filtered_1 <- converted |> 
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


write_csv(normalized, snakemake@output[[1]])