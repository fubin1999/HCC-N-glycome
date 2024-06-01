library(tidyverse)
library(rstatix)

source("src/utils/emeans_posthoc.R")

# residues <- read_csv("results/data/residues/glycan_residues.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv") %>%
#   select(sample, age, sex)

residues <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clinical <- read_csv(snakemake@input[[3]]) %>%
  select(sample, age, sex)

data <- residues %>%
  pivot_longer(-sample, names_to = "residue", values_to = "value") %>%
  inner_join(groups, by = "sample") %>%
  inner_join(clinical, by = "sample") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

ancova_result <- data %>%
  group_by(residue) %>%
  anova_test(value ~ group + age + sex, white.adjust = TRUE) %>%
  adjust_pvalue() %>%
  add_significance("p.adj") %>%
  as_tibble() %>%
  select(-`p<.05`)

post_hoc_result <- data %>%
  post_hoc(value ~ group + age + sex, group, residue)

write_csv(ancova_result, snakemake@output[[1]])
write_csv(post_hoc_result, snakemake@output[[2]])