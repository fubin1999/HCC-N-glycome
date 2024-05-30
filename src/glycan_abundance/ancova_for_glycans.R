source("renv/activate.R")

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
  inner_join(clinical |> select(sample, sex, age), by = "sample") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  mutate(log_value = log2(value))

# ANCOVA-----
ancova_result <- data |> 
  group_by(glycan) |>
  anova_test(
    log_value ~ group + age + sex,
    white.adjust = TRUE,
    effect.size = "pes"
  ) |>
  adjust_pvalue(method = "BH") |> 
  as_tibble()

# Post-hoc-----
diff_glycans <- ancova_result |> 
  filter(p.adj < 0.05, Effect == "group") |>
  pull(glycan)

posthoc_result <- data %>%
  filter(glycan %in% diff_glycans) %>%
  post_hoc(log_value ~ group + sex + age, group, glycan) %>%
  add_significance(p.col = "p.adj")

# Save results-----
write_csv(ancova_result, snakemake@output[["ancova"]])
write_csv(posthoc_result, snakemake@output[["posthoc"]])