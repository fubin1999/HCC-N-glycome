library(tidyverse)
library(rstatix)
library(patchwork)

trait_data <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clinical <- read_csv(snakemake@input[[3]]) %>%
  select(-age, -sex)

data <- trait_data %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") %>%
  inner_join(groups, by = "sample") %>%
  inner_join(
    pivot_longer(clinical, -sample, names_to = "clinical", values_to = "clinical_value"),
    by = "sample", relationship = "many-to-many"
  )

calcu_corr <- function (data) {
  data %>%
    group_by(trait, clinical) %>%
    cor_test(value, clinical_value, method = "spearman") %>%
    adjust_pvalue() %>%
    as_tibble() %>%
    select(trait, clinical, cor, p, p.adj)
}

full_cor_result <- calcu_corr(data)

sub_cor_result <- data %>%
  nest_by(group) %>%
  mutate(cor_result = list(calcu_corr(data))) %>%
  select(-data) %>%
  unnest(cor_result)

write_csv(full_cor_result, snakemake@output[[1]])
write_csv(sub_cor_result, snakemake@output[[2]])
