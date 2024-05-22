source("renv/activate.R")

library(tidyverse)

# Read data-----
# glycan_data <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
glycan_data <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

data <- glycan_data %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC")

# Calculate Fold Change-----
FC_result <- data %>%
  summarise(mean_value = mean(value), .by = c(glycan, group)) %>%
  pivot_wider(names_from = group, values_from = mean_value) %>%
  mutate(HCC_CHB = HCC / CHB, HCC_HC = HCC / HC, HCC_LC = HCC / LC) %>%
  select(glycan, HCC_CHB, HCC_HC, HCC_LC) %>%
  pivot_longer(-glycan, names_to = "comparison", values_to = "FC") %>%
  separate_wider_delim(comparison, delim = "_", names = c("group1", "group2"))

# Save result-----
write_csv(FC_result, snakemake@output[[1]])