library(tidyverse)

# Load data-----
abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

# Total Abundance Normalization-----
selected <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "abundance") %>%
  summarise(missing_prop = mean(is.na(abundance)), .by = glycan) %>%
  filter(missing_prop < 0.01) %>%
  pull(glycan)

normed <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "abundance") %>%
  filter(glycan %in% selected) %>%
  mutate(abundance = abundance / sum(abundance, na.rm = TRUE) * 100, .by = sample) %>%
  pivot_wider(names_from = glycan, values_from = abundance)

prepared <- groups %>%
  inner_join(normed, by = "sample") %>%
  filter(group != "QC")

# Train-Test Split-----
set.seed(123)
train_data <- prepared %>%
  slice_sample(prop = 0.7, by = group)
test_data <- prepared %>%
  filter(!sample %in% train_data$sample)

# Save Data-----
write_csv(train_data, snakemake@output[[1]])
write_csv(test_data, snakemake@output[[2]])
