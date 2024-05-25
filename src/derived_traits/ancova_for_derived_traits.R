source("renv/activate.R")

library(tidyverse)
library(rstatix)

# Read data-----
# trait_data <- read_csv("results/data/derived_traits/derived_traits.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv") %>%
#   select(sample, sex, age)

trait_data <- read_csv(snakemake@input[["traits"]])
groups <- read_csv(snakemake@input[["groups"]])
clinical <- read_csv(snakemake@input[["clinical"]]) %>%
  select(sample, sex, age)

data <- trait_data %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") %>%
  inner_join(groups, by = "sample") %>%
  inner_join(clinical, by = "sample") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

# Perform ANCOVA-----
ancova_result <- data |>
  group_by(trait) |>
  anova_test(
    value ~ group + age + sex,
    white.adjust = TRUE,
    effect.size = "pes"
  ) |>
  adjust_pvalue(method = "BH") |>
  as_tibble()

# Post-hoc-----
diff_traits <- ancova_result |>
  filter(p.adj < 0.05, Effect == "group") |>
  pull(trait)

posthoc_result <- data %>%
  filter(trait %in% diff_traits) %>%
  group_by(trait) %>%
  nest() %>%
  mutate(
    lm_model = map(data, ~ lm(value ~ group + sex + age, data = .x)),
    emm = map(lm_model, ~ emmeans::emmeans(.x, ~ group)),
    pairwise_comparison = map(emm, ~ emmeans::contrast(.x, "pairwise", adjust = "tukey")),
    result_df = map(pairwise_comparison, ~ as_tibble(.x))
  ) %>%
  select(trait, result_df) %>%
  unnest(cols = result_df) %>%
  separate_wider_delim(contrast, delim = " - ", names = c("group1", "group2")) %>%
  rename(p.adj = p.value) %>%
  add_significance(p.col = "p.adj")

# Save results-----
write_csv(ancova_result, snakemake@output[[1]])
write_csv(posthoc_result, snakemake@output[[2]])