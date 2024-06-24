library(tidyverse)
library(rstatix)


# traits <- read_csv("results/data/prepared/filtered_derived_traits.csv")
# subtypes <- read_csv("results/data/subtype/cc_result.csv")

traits <- read_csv(snakemake@input[[1]])
subtypes <- read_csv(snakemake@input[[2]])

data <- traits %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") %>%
  right_join(subtypes, by = "sample") %>%
  mutate(class = factor(class))

t_test_result <- data %>%
  group_by(trait) %>%
  t_test(value ~ class) %>%
  adjust_pvalue(method = "BH") %>%
  as_tibble() %>%
  select(trait, statistic, p, p.adj)

eff_size <- data %>%
  group_by(trait) %>%
  cohens_d(value ~ class) %>%
  select(trait, cohens_d = effsize, magnitude)

final_result <- t_test_result %>%
  left_join(eff_size, by = "trait")

write_csv(final_result, snakemake@output[[1]])