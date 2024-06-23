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

HCC_data <- data %>%
  filter(group == "HCC")

calcu_corr <- function (data) {
  data %>%
    group_by(trait, clinical) %>%
    cor_test(value, clinical_value, method = "spearman") %>%
    adjust_pvalue() %>%
    as_tibble() %>%
    select(trait, clinical, cor, p, p.adj)
}

full_cor_result <- calcu_corr(data)
HCC_cor_result <- calcu_corr(HCC_data)

write_csv(full_cor_result, snakemake@output[[1]])
write_csv(HCC_cor_result, snakemake@output[[2]])
