library(tidyverse)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# mp_table <- read_csv("results/data/derived_traits/meta_properties.csv")
abundance <- read_csv(snakemake@input[[1]])
mp_table <- read_csv(snakemake@input[[2]])

residue_table <- mp_table %>%
  select(glycan, nF, nS, nM, nG, nN) %>%
  pivot_longer(-glycan, names_to = "residue", values_to = "n", names_prefix = "n") %>%
  mutate(
    n = as.integer(n),
    residue = as.factor(residue)
  )

total_abundance <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  summarise(total_abundance = sum(value), .by = "sample")

result <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  left_join(residue_table, by = "glycan", relationship = "many-to-many") %>%
  summarise(value = sum(value * n), .by = c(sample, residue)) %>%
  left_join(total_abundance, by = "sample") %>%
  mutate(value = value / total_abundance) %>%
  select(-total_abundance) %>%
  pivot_wider(names_from = "residue", values_from = "value")

write_csv(result, snakemake@output[[1]])