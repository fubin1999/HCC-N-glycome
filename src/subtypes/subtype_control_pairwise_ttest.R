library(tidyverse)
library(rstatix)

# Read and prepare data-----
# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# subtypes <- read_csv("results/data/subtypes/consensus_cluster_result.csv")

abundance <- read_csv(snakemake@input[["abundance"]])
groups <- read_csv(snakemake@input[["groups"]])
subtypes <- read_csv(snakemake@input[["subtypes"]])

prepared <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  left_join(subtypes %>% mutate(subtype = paste0("S", class), .keep = "unused"), by = "sample") %>%
  mutate(group = if_else(group == "HCC", subtype, group)) %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "S1", "S2", "S3"))) %>%
  select(-subtype)

# Pairwise t-test-----
ttest_result <- prepared %>%
  mutate(value = log(value)) %>%
  group_by(glycan) %>%
  t_test(
    value ~ group,
    comparisons = list(
      c("S1", "HC"), c("S1", "CHB"), c("S1", "LC"),
      c("S2", "HC"), c("S2", "CHB"), c("S2", "LC"),
      c("S3", "HC"), c("S3", "CHB"), c("S3", "LC")
    ),
    p.adjust.method = "none",
  ) %>%
  select(-c(p.adj, p.adj.signif)) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  as_tibble() %>%
  select(-c(.y., n1, n2, df, p))

# Pairwise fold change-----
grouped_means <- prepared %>%
  summarise(mean_value = mean(value), .by = c("glycan", "group"))

fold_change_result <- expand_grid(
  glycan = setdiff(colnames(abundance), "sample"),
  group1 = c("HC", "CHB", "LC"),
  group2 = c("S1", "S2", "S3")
) %>%
  left_join(grouped_means, by = c("glycan" = "glycan", "group1" = "group")) %>%
  rename(mean_group_1 = mean_value) %>%
  left_join(grouped_means, by = c("glycan" = "glycan", "group2" = "group")) %>%
  rename(mean_group_2 = mean_value) %>%
  mutate(fold_change = mean_group_2 / mean_group_1) %>%
  select(-c(mean_group_1, mean_group_2))

# Combine and save results-----
result <- ttest_result %>%
  left_join(fold_change_result, by = c("glycan", "group1", "group2"))

write_csv(result, snakemake@output[[1]])