library(tidyverse)

raw_data <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clinical <- read_csv(snakemake@input[[3]])
plates <- read_csv(snakemake@input[[4]])
chosen_samples <- read_csv(snakemake@input[[5]])$sample

# raw_data <- read_csv("results/data/prepared/raw_abundance.csv")
# groups <- read_csv("results/data/prepared/unfiltered_groups.csv")
# clinical <- read_csv("results/data/prepared/unfiltered_clinical.csv")
# plates <- read_csv("data/cohort_GD1/plates.csv")
# chosen_samples <- read_csv("data/cohort_GD1/chosen_samples.csv")$sample

# 1. Convert to long-----
long_data <- raw_data |>
  filter(sample %in% chosen_samples) |>
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
glycans_to_keep <- filtered_1 %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  group_by(glycan, group) %>%
  summarise(detect_rate = mean(!is.na(value)), .groups = "drop") %>%
  filter(detect_rate > 0.5) %>%
  distinct(glycan) %>%
  pull(glycan)

filtered_2 <- filtered_1 %>%
  filter(glycan %in% glycans_to_keep)

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

# 6. Batch effect correction-----
library(sva)

to_be_combated <- normalized |>
  pivot_wider(names_from = sample, values_from = value) %>%
  column_to_rownames(var = "glycan") %>%
  as.matrix() %>%
  log()
batch_data <- plates %>%
  mutate(group = str_extract(raw_sample, "^[A-Z]+")) %>%
  column_to_rownames("sample") %>%
  .[colnames(to_be_combated), c("plate", "group")] %>%
  mutate(plate = as.factor(plate), group = as.factor(group))
mod_combat <- model.matrix(~group, data = batch_data)
combat_data <- ComBat(dat = to_be_combated, batch = batch_data$plate, mod = mod_combat, ref.batch = 1)
final_data <- combat_data %>%
  exp() %>%
  t() %>% as.data.frame() %>% rownames_to_column("sample") %>%
  as_tibble()

# Filter groups and clinical-----
filtered_groups <- groups %>%
  semi_join(final_data, by = "sample")
filtered_clinical <- clinical %>%
  semi_join(filtered_groups %>% filter(group != "QC"), by = "sample")

# Save result-----
write_csv(final_data, snakemake@output[[1]])
write_csv(filtered_groups, snakemake@output[[2]])
write_csv(filtered_clinical, snakemake@output[[3]])