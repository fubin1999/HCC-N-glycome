library(tidyverse)
library(rstatix)
library(patchwork)

# glycan_data <- read_csv("results/data/prepared/processed_abundance.csv")
# trait_data <- read_csv("results/data/prepared/filtered_derived_traits.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv") %>%
#   select(-age, -sex)

glycan_data <- read_csv(snakemake@input[["glycans"]])
trait_data <- read_csv(snakemake@input[["traits"]])
groups <- read_csv(snakemake@input[["groups"]])
clinical <- read_csv(snakemake@input[["clinical"]]) %>%
  select(-age, -sex)

prepare_data <- function (data) {
  data %>%
    pivot_longer(-sample, names_to = "feature", values_to = "value") %>%
    inner_join(groups, by = "sample") %>%
    inner_join(
      pivot_longer(clinical, -sample, names_to = "clinical", values_to = "clinical_value"),
      by = "sample", relationship = "many-to-many"
    )
}

glycan_data_prepared <- prepare_data(glycan_data)
trait_data_prepared <- prepare_data(trait_data)

calcu_corr <- function (data) {
  data %>%
    group_by(feature, clinical) %>%
    cor_test(value, clinical_value, method = "spearman") %>%
    adjust_pvalue() %>%
    as_tibble() %>%
    select(feature, clinical, cor, p, p.adj)
}

glycan_full_cor_result <- calcu_corr(glycan_data_prepared) %>%
  rename(glycan = feature)
trait_full_cor_result <- calcu_corr(trait_data_prepared) %>%
  rename(trait = feature)

glycan_sub_cor_result <- glycan_data_prepared %>%
  nest_by(group) %>%
  mutate(cor_result = list(calcu_corr(data))) %>%
  select(-data) %>%
  unnest(cor_result) %>%
  rename(glycan = feature)

trait_sub_cor_result <- trait_data_prepared %>%
  nest_by(group) %>%
  mutate(cor_result = list(calcu_corr(data))) %>%
  select(-data) %>%
  unnest(cor_result) %>%
  rename(trait = feature)

write_csv(glycan_full_cor_result, snakemake@output[[1]])
write_csv(glycan_sub_cor_result, snakemake@output[[2]])
write_csv(trait_full_cor_result, snakemake@output[[3]])
write_csv(trait_sub_cor_result, snakemake@output[[4]])