library(tidyverse)
library(rstatix)

# raw_abundance <- read_csv("results/data/prepared/raw_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
raw_abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

detection <- raw_abundance %>%
  right_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  pivot_longer(-c(sample, group), names_to = "glycan", values_to = "value") %>%
  mutate(detected = !is.na(value)) %>%
  group_by(group, glycan, detected) %>%
  count() %>%
  ungroup() %>%
  complete(group, glycan, detected, fill = list(n = 0))

tidy_fisher_test <- function (data) {
  data %>%
    pivot_wider(names_from = group, values_from = n) %>%
    column_to_rownames("detected") %>%
    as.matrix() %>% as.table() %>%
    pairwise_fisher_test(detailed = TRUE, p.adjust.method = "BH") %>%
    as_tibble() %>%
    select(group1, group2, estimate, conf.low, conf.high, p)
}

fisher_result <- detection %>%
  nest_by(glycan) %>%
  mutate(fisher_result = list(tidy_fisher_test(data))) %>%
  select(-data) %>%
  unnest(fisher_result) %>%
  adjust_pvalue()

write_csv(fisher_result, snakemake@output[[1]])