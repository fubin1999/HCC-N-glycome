source("renv/activate.R")

library(tidyverse)
library(rstatix)

# Read data-----
glycan_abundance <- read_csv(snakemake@input[["abundance"]])
groups <- read_csv(snakemake@input[["groups"]])
clinical <- read_csv(snakemake@input[["clinical"]])

data <- glycan_abundance %>%
  inner_join(groups, by = "sample") %>%
  inner_join(clinical |> select(sample, sex, age), by = "sample")

# ANCOVA-----
ancova_result <- data |> 
  mutate(log_value = log2(value)) |> 
  group_by(glycan) |>
  anova_test(log_value ~ group + age + sex, white.adjust = TRUE) |> 
  adjust_pvalue(method = "BH") |> 
  as_tibble() |> 
  filter(Effect == "group") |> 
  select(-Effect)

# Post-hoc-----
diff_glycans <- ancova_result |> 
  filter(p.adj < 0.05) |> 
  pull(glycan)

posthoc_result <- data |>
  filter(glycan %in% diff_glycans) |>
  mutate(log_value = log2(value)) |>
  group_by(glycan) |>
  tukey_hsd(log_value ~ group, p.adjust.method = "BH") |>
  as_tibble() |> 
  filter(term == "group") |>
  select(-term)

# Save results-----
write_csv(ancova_result, snakemake@output[["ancova"]])
write_csv(posthoc_result, snakemake@output[["posthoc"]])