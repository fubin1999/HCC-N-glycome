library(tidyverse)
library(rstatix)

source("src/utils/emeans_posthoc.R")

# Read data-----
motifs <- read_csv("results/data/motifs/motifs.csv")
groups <- read_csv("results/data/prepared/groups.csv")
clinical <- read_csv("results/data/prepared/clinical.csv")

motifs <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clinical <- read_csv(snakemake@input[[3]])

data <- motifs %>%
  pivot_longer(-sample, names_to = "motif", values_to = "value") %>%
  inner_join(groups, by = "sample") %>%
  inner_join(clinical |> select(sample, sex, age), by = "sample") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  mutate(log_value = log2(value * 100 + 1))

low_variance_motifs <- data %>%
  summarise(var = sd(log_value), .by = motif) %>%
  filter(var == 0) %>%
  pull(motif)

data <- data %>%
  filter(!(motif %in% low_variance_motifs))

# ANCOVA-----
ancova_result <- data |>
  group_by(motif) |>
  anova_test(
    log_value ~ group + age + sex,
    white.adjust = TRUE,
    effect.size = "pes"
  ) |>
  adjust_pvalue(method = "BH") |>
  as_tibble()

# Post-hoc-----
diff_motifs <- ancova_result |>
  filter(p.adj < 0.05, Effect == "group") |>
  pull(motif)

posthoc_result <- data %>%
  filter(motif %in% diff_motifs) %>%
  post_hoc(log_value ~ group + sex + age, group, motif) %>%
  add_significance(p.col = "p.adj")

# Save results-----
write_csv(ancova_result, snakemake@output[[1]])
write_csv(posthoc_result, snakemake@output[[2]])