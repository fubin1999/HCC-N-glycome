source("renv/activate.R")

library(tidyverse)
library(rstatix)

source("src/utils/emeans_posthoc.R")

# Read data-----
trait_data <- read_csv("results/data/derived_traits/filtered_derived_traits.csv")
groups <- read_csv("results/data/prepared/groups.csv")
clinical <- read_csv("results/data/prepared/clinical.csv") %>%
  select(sample, sex, age)

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
  post_hoc(value ~ group + sex + age, group, trait) %>%
  add_significance(p.col = "p.adj")

# Save results-----
write_csv(ancova_result, snakemake@output[[1]])
write_csv(posthoc_result, snakemake@output[[2]])