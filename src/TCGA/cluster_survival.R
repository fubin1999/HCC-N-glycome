library(SummarizedExperiment)
library(survival)
library(survminer)
library(tidyverse)
library(patchwork)

# load("results/data/TCGA/prepared_data.rda")
# cc_result <- read_csv("results/data/TCGA/consensus_cluster_result.csv")

load(snakemake@input[[1]])
cc_result <- read_csv(snakemake@input[[2]])

sample_info <- as_tibble(colData(data))
HCC_sample_info <- filter(sample_info, shortLetterCode == "TP")

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
  left_join(cc_result, by = "barcode")

survival_result <- survfit(Surv(os, vital_status) ~ class, data = survival_data)
p_value <- surv_pvalue(survival_result, data = survival_data, method = "sqrtN")

cox_model <- coxph(Surv(os, vital_status) ~ class, data = survival_data)

plots <- ggsurvplot(
  survival_result,
  data = survival_data,
  legend.labs = c("Cluster1", "Cluster2"),
  palette = "npg",
  risk.table = TRUE,
  surv.median.line = "hv",
  pval = TRUE,
  pval.method = TRUE,
  log.rank.weights = "sqrtN",
  pval.coord = c(2100, 0.85),
  pval.method.coord = c(2100, 0.95),
)

final_p <- plots$plot / plots$table + plot_layout(heights = c(4, 1))
# tgutil::ggpreview(plot = final_p, width = 5, height = 7)
ggsave(snakemake@output[[1]], plot = final_p, width = 5, height = 7)