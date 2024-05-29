source("renv/activate.R")

library(tidyverse)
library(gt)
library(webshot2)

# metrics_file <- "results/data/ml/model_performance.csv"
metrics_file <- snakemake@input[[1]]
metrics <- read_csv(metrics_file) |> 
  mutate(
    comparison = factor(comparison, levels = rev(c("Control/HCC", "HC/HCC", "CHB/HCC", "LC/HCC"))),
    metric = case_match(
      metric,
      "accuracy" ~ "Accuracy",
      "sensitivity" ~ "Sensitivity",
      "specificity" ~ "Specificity",
      "f1_score" ~ "F1 Score",
      "roc_auc" ~ "ROC AUC",
      "pr_auc" ~ "PR AUC",
    ),
    metric = factor(metric, levels = c("Accuracy", "Sensitivity", "Specificity", "F1 Score", "ROC AUC", "PR AUC"))
  ) |> 
  pivot_wider(names_from = metric, values_from = score)

table <- gt(metrics, rowname_col = "comparison") |> 
  tab_header(title = md("**Model Performance**")) |>
  fmt_number(columns = 2:7, decimals = 3) |> 
  tab_footnote("ROC AUC: Receiver Operating Characteristic Area Under the Curve") |> 
  tab_footnote("PR AUC: Precision-Recall Area Under the Curve")
# gtsave(table, "results/figures/ml/model_performance.pdf")
gtsave(table, snakemake@output[[1]])
