library(tidyverse)
library(rstatix)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# clusters <- read_csv("results/data/subtypes/consensus_cluster_result.csv")

abundance <- read_csv(snakemake@input[[1]])
clusters <- read_csv(snakemake@input[[2]])

data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  mutate(value = log(value)) %>%
  right_join(clusters, by = "sample") %>%
  mutate(class = as.factor(class))

anova_result <- data %>%
  group_by(glycan) %>%
  anova_test(value ~ class) %>%
  adjust_pvalue(method = "BH") %>%
  as_tibble()

post_hoc_result <- data %>%
  group_by(glycan) %>%
  tukey_hsd(value ~ class) %>%
  as_tibble()

write_csv(anova_result, snakemake@output[[1]])
write_csv(post_hoc_result, snakemake@output[[2]])