library(tidyverse)
library(rstatix)
library(broom)

source("src/utils/emeans_posthoc.R")

# Read data-----
# glycan_abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")

glycan_abundance <- read_csv(snakemake@input[["abundance"]])
groups <- read_csv(snakemake@input[["groups"]])
clinical <- read_csv(snakemake@input[["clinical"]])

data <- glycan_abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  inner_join(groups, by = "sample") %>%
  inner_join(clinical |> select(sample, sex, age, ALBI_score, AAR), by = "sample") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  mutate(log_value = log2(value))

# ANCOVA-----
ancova_result <- data |> 
  group_by(glycan) |>
  anova_test(
    log_value ~ age * sex + group,
    white.adjust = TRUE,
    effect.size = "pes"
  ) |>
  adjust_pvalue(method = "BH") |> 
  as_tibble()

ancova_result_lf_adjusted <- data |>
  group_by(glycan) |>
  anova_test(
    log_value ~ age * sex + ALBI_score + group,
    white.adjust = TRUE,
    effect.size = "pes"
  ) |>
  adjust_pvalue(method = "BH") |>
  as_tibble()

# Post-hoc-----
posthoc_result <- data %>%
  post_hoc(log_value ~ age * sex + group, group, glycan) %>%
  add_significance(p.col = "p.adj")

# Post-hoc-----
posthoc_result_lf_adjusted <- data %>%
  post_hoc(log_value ~ age * sex + ALBI_score + group, group, glycan) %>%
  add_significance(p.col = "p.adj")

# Save results-----
write_csv(ancova_result, snakemake@output[["ancova"]])
write_csv(posthoc_result, snakemake@output[["posthoc"]])
write_csv(ancova_result_lf_adjusted, snakemake@output[["ancova_lf_adjusted"]])
write_csv(posthoc_result_lf_adjusted, snakemake@output[["posthoc_lf_adjusted"]])