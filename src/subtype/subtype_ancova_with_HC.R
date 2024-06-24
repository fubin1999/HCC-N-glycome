library(tidyverse)
library(rstatix)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")
# subtypes <- read_csv("results/data/subtype/cc_result.csv")

abundance <- read_csv(snakemake@input[["abundance"]])
groups <- read_csv(snakemake@input[["groups"]])
clinical <- read_csv(snakemake@input[["clinical"]])
subtypes <- read_csv(snakemake@input[["subtypes"]])

data <- groups %>%
  left_join(subtypes, by = "sample") %>%
  left_join(clinical %>% select(sample, sex, age), by = "sample") %>%
  left_join(abundance, by = "sample") %>%
  mutate(group = case_when(
    group == "HC" ~ "HC",
    class == 1 ~ "HCC_S1",
    class == 2 ~ "HCC_S2"
  )) %>%
  filter(!is.na(group)) %>%
  select(-class) %>%
  pivot_longer(starts_with("H"), names_to = "glycan", values_to = "value")

S1_result <- data %>%
  filter(group != "HCC_S2") %>%
  mutate(value = log(value + 1)) %>%
  group_by(glycan) %>%
  anova_test(value ~ sex * age + group, white.adjust = TRUE) %>%
  as_tibble()

S2_result <- data %>%
  filter(group != "HCC_S1") %>%
  mutate(value = log(value + 1)) %>%
  group_by(glycan) %>%
  anova_test(value ~ sex * age + group, white.adjust = TRUE) %>%
  as_tibble()

ancova_result <- bind_rows(list(HCC_S1 = S1_result, HCC_S2 = S2_result), .id = "group2") %>%
  mutate(group1 = "HC") %>%
  filter(Effect == "group") %>%
  adjust_pvalue(method = "BH") %>%
  select(glycan, group1, group2, F, p, p.adj)

S1_fold_change <- data %>%
  filter(group != "HCC_S2") %>%
  group_by(glycan, group) %>%
  summarise(mean = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = "group", values_from = "mean") %>%
  mutate(fold_change = HCC_S1 / HC) %>%
  select(glycan, fold_change)

S2_fold_change <- data %>%
  filter(group != "HCC_S1") %>%
  group_by(glycan, group) %>%
  summarise(mean = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = "group", values_from = "mean") %>%
  mutate(fold_change = HCC_S2 / HC) %>%
  select(glycan, fold_change)

fold_change <- bind_rows(list(HCC_S1 = S1_fold_change, HCC_S2 = S2_fold_change), .id = "group2") %>%
  mutate(group1 = "HC") %>%
  select(glycan, group1, group2, fold_change)

final_result <- ancova_result %>%
  left_join(fold_change, by = c("glycan", "group1", "group2"))

write_csv(final_result, snakemake@output[[1]])