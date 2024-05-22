source("renv/activate.R")

library(tidyverse)
library(rstatix)

# Read data-----
derived_traits <- read_csv(snakemake@input[["derived_traits"]])
groups <- read_csv(snakemake@input[["groups"]])
clinical <- read_csv(snakemake@input[["clinical"]])

data <- derived_traits %>%
  inner_join(groups, by = "sample") %>%
  inner_join(clinical |> select(sample, sex, age), by = "sample") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

# ANCOVA-----
ancova_result <- data |> 
  group_by(trait) |>
  anova_test(value ~ group + age + sex, white.adjust = TRUE) |> 
  adjust_pvalue(method = "BH") |> 
  as_tibble() |> 
  filter(Effect == "group") |> 
  select(-Effect)

# Post-hoc-----
diff_traits <- ancova_result |> 
  filter(p.adj < 0.05) |> 
  pull(trait)

posthoc_result <- data |>
  filter(trait %in% diff_traits) |>
  group_by(trait) |>
  tukey_hsd(value ~ group, p.adjust.method = "BH") |>
  as_tibble() |> 
  filter(term == "group") |>
  select(-term)

# Save results-----
write_csv(ancova_result, snakemake@output[["ancova"]])
write_csv(posthoc_result, snakemake@output[["posthoc"]])