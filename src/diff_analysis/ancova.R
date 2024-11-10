library(tidyverse)
library(rstatix)
library(broom)

source("src/utils/emeans_posthoc.R")

# Read data-----
# data <- read_csv("results/data/prepared/filtered_derived_traits.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")
# var_name <- "trait"
# cov_albi <- FALSE

data <- read_csv(snakemake@input[["data"]])
groups <- read_csv(snakemake@input[["groups"]])
clinical <- read_csv(snakemake@input[["clinical"]])
var_name <- snakemake@params[["var_name"]]
cov_albi <- snakemake@params[["cov_albi"]]

data <- data %>%
  pivot_longer(-sample, names_to = var_name, values_to = "value") %>%
  inner_join(groups, by = "sample") %>%
  inner_join(clinical |> select(sample, sex, age, ALBI_score), by = "sample") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

if (var_name == "glycan") {
  data <- data %>% mutate(value = log2(value))
}

if (cov_albi) {
  formula <- value ~ age * sex + ALBI_score + group
} else {
  formula <- value ~ age * sex + group
}

# ANCOVA-----
ancova_result <- data |> 
  group_by(.data[[var_name]]) |>
  anova_test(formula, white.adjust = TRUE, effect.size = "pes") |>
  adjust_pvalue(method = "BH") |> 
  as_tibble()

# Post-hoc-----
posthoc_result <- data %>%
  post_hoc(formula, group, !!rlang::sym(var_name)) %>%
  add_significance(p.col = "p.adj")

# Save results-----
write_csv(ancova_result, snakemake@output[[1]])
write_csv(posthoc_result, snakemake@output[[2]])
