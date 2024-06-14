library(SummarizedExperiment)
library(survival)
library(survminer)
library(tidyverse)
library(broom)
library(rstatix)


# load("results/data/TCGA/prepared_data.rda")
# glycogenes <- read_csv("data/glycogenes.csv")$gene_name

load(snakemake@input[[1]])
glycogenes <- read_csv(snakemake@input[[2]])$gene_name

sample_info <- as_tibble(colData(data))
gene_info <- as_tibble(rowData(data))
glyco_gene_info <- gene_info %>% filter(gene_name %in% glycogenes)
expr_mat <- assay(data, "tpm_unstrand")

HCC_sample_info <- filter(sample_info, shortLetterCode == "TP")
HCC_glycogene_expr_mat <- expr_mat[
  glyco_gene_info$gene_id,
  HCC_sample_info$barcode
]

tidy_expr_data <- HCC_glycogene_expr_mat %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  pivot_longer(-gene_id, names_to = "barcode", values_to = "expr") %>%
  left_join(glyco_gene_info %>% select(gene_name, gene_id), by = "gene_id") %>%
  select(-gene_id)

survival_data <- HCC_sample_info %>%
  select(barcode, days_to_last_follow_up, days_to_death, vital_status, gender) %>%
  mutate(gender = factor(gender)) %>%
  mutate(vital_status = if_else(vital_status == "Dead", 1, 0)) %>%
  mutate(os = if_else(
    vital_status == 1,
    days_to_death,
    days_to_last_follow_up
  )) %>%
  select(barcode, os, vital_status, gender) %>%
  filter(!is.na(os)) %>%
  left_join(tidy_expr_data, by = "barcode")

survival_results <- survival_data %>%
  nest_by(gene_name) %>%
  mutate(
    fit = list(coxph(Surv(os, vital_status) ~ gender + expr, data = data)),
    fit_stats = list(tidy(fit, conf.int = TRUE, exponentiate = TRUE))
  ) %>%
  select(gene_name, fit_stats) %>%
  unnest(fit_stats) %>%
  filter(term == "expr") %>%
  select(-term) %>%
  adjust_pvalue(p.col = "p.value", method = "BH")

write_csv(survival_results, snakemake@output[[1]])