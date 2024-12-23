library(tidyverse)
library(pROC)


# Read data-----
# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])


# Prepare data-----
with_groups <- abundance %>% 
  pivot_longer(-sample, names_to = "glycan", values_to = "abundance") %>% 
  left_join(groups, by = "sample") %>% 
  filter(!group == "QC")

HCC_vs_rest <- with_groups %>% 
  mutate(group = if_else(group == "HCC", "HCC", "Control"))
HCC_vs_HC <- with_groups %>% 
  filter(group %in% c("HCC", "HC"))
HCC_vs_CHB <- with_groups %>% 
  filter(group %in% c("HCC", "CHB"))
HCC_vs_LC <- with_groups %>%
  filter(group %in% c("HCC", "LC"))
prepared <- bind_rows(list(
  HCC_vs_rest = HCC_vs_rest,
  HCC_vs_HC = HCC_vs_HC,
  HCC_vs_CHB = HCC_vs_CHB,
  HCC_vs_LC = HCC_vs_LC
), .id = "comparison")


# Calculate AUC-----
roc_aucs <- prepared %>% 
  nest_by(comparison, glycan) %>% 
  mutate(
    roc = list(roc(group ~ abundance, data = data, ci = TRUE, ci.method = "delong")),
    auc = as.numeric(auc(roc)),
    ci_lower = as.numeric(ci(roc)[[1]]),
    ci_upper = as.numeric(ci(roc)[[3]])
  ) %>% 
  select(comparison, glycan, auc, ci_lower, ci_upper)

# Save results-----
write_csv(roc_aucs, snakemake@output[[1]])