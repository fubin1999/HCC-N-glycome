source("renv/activate.R")

library(tidyverse)

# Read data-----
# trait_data <- read_csv("results/data/derived_traits/filtered_derived_traits.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
trait_data <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

data <- trait_data %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") |>
    left_join(groups, by = "sample") %>%
  filter(group != "QC")

# Calculate Fold Change-----
FC_result <- data %>%
  summarise(mean_value = mean(value), .by = c(trait, group)) %>%
  pivot_wider(names_from = group, values_from = mean_value) %>%
  mutate(
    CHB_HCC = HCC / CHB, HC_HCC = HCC / HC, LC_HCC = HCC / LC,
    HC_LC = LC / HC, CHB_LC = LC / CHB, HC_CHB = CHB / HC
  ) %>%
  select(trait, CHB_HCC, HC_HCC, LC_HCC, HC_LC, CHB_LC, HC_CHB) %>%
  pivot_longer(-trait, names_to = "comparison", values_to = "FC") %>%
  separate_wider_delim(comparison, delim = "_", names = c("group1", "group2"))

# Save result-----
write_csv(FC_result, snakemake@output[[1]])