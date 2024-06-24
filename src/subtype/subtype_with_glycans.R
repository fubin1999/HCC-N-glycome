library(tidyverse)
library(rstatix)


# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# subtypes <- read_csv("results/data/subtype/cc_result.csv")

abundance <- read_csv(snakemake@input[[1]])
subtypes <- read_csv(snakemake@input[[2]])

data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  right_join(subtypes, by = "sample") %>%
  mutate(class = factor(class))

t_test_result <- data %>%
  mutate(value = log(value + 1)) %>%
  group_by(glycan) %>%
  t_test(value ~ class) %>%
  adjust_pvalue(method = "BH") %>%
  as_tibble() %>%
  select(glycan, statistic, p, p.adj)

fold_change <- data %>%
  group_by(glycan, class) %>%
  summarise(mean = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = class, values_from = mean, names_prefix = "S") %>%
  mutate(fold_change = S1 / S2) %>%
  select(glycan, fold_change)

eff_size <- data %>%
  mutate(value = log(value + 1)) %>%
  group_by(glycan) %>%
  cohens_d(value ~ class) %>%
  as_tibble() %>%
  select(glycan, cohens_d = effsize, magnitude)

final_result <- t_test_result %>%
  left_join(fold_change, by = "glycan") %>%
  left_join(eff_size, by = "glycan")

write_csv(final_result, snakemake@output[[1]])