library(SummarizedExperiment)
library(survival)
library(survminer)
library(tidyverse)
library(patchwork)

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
  select(barcode, days_to_last_follow_up, days_to_death, vital_status) %>%
  mutate(vital_status = if_else(vital_status == "Dead", 1, 0)) %>%
  mutate(os = if_else(
    vital_status == 1,
    days_to_death,
    days_to_last_follow_up
  )) %>%
  select(barcode, os, vital_status) %>%
  filter(!is.na(os)) %>%
  left_join(tidy_expr_data, by = "barcode") %>%
  pivot_wider(names_from = gene_name, values_from = expr) %>%
  surv_cutpoint(
    time = "os", event = "vital_status",
    variables = glyco_gene_info$gene_name
  ) %>%
  surv_categorize() %>%
  as_tibble()

survival_results <- survival_data %>%
  pivot_longer(-c(os, vital_status), names_to = "gene_name", values_to = "group") %>%
  nest_by(gene_name) %>%
  mutate(fit = list(survfit(Surv(os, vital_status) ~ group, data = data)))

p_values <- survival_results %>%
  mutate(p_value = surv_pvalue(fit, data)) %>%
  select(gene_name, p_value) %>%
  unnest(p_value)

sig_genes <- p_values %>%
  filter(pval < 0.05) %>%
  pull(gene_name)

surv_plots <- survival_results %>%
  filter(gene_name %in% sig_genes) %>%
  mutate(plot = list(ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    conf.int = TRUE,
    legend.labs = c("High", "Low"),
    palette = "npg",
    risk.table = TRUE,
    surv.median.line = "hv",
    conf.int.alpha = 0.1,
    pval.coord = c(2400, 0.95)
  ))) %>%
  select(gene_name, plot)

dirpath <- snakemake@output[[1]]
walk2(
  str_glue("{dirpath}/{surv_plots$gene_name}.pdf"),
  surv_plots$plot,
  ~ ggsave(.x, survminer:::.build_ggsurvplot(.y), create.dir = TRUE)
)