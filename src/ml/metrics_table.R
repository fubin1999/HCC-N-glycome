library(tidyverse)
library(gt)
library(webshot2)

# metrics_file <- "results/data/ml/model_performance.csv"
metrics_file <- snakemake@input[[1]]

metrics <- read_csv(metrics_file) |> 
  mutate(
    comparison = str_replace(comparison, "/", " / "),
    comparison = factor(comparison, levels = c("HC / HCC", "CHB / HCC", "LC / HCC", "Control / HCC")),
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
  pivot_wider(names_from = metric, values_from = score) %>%
  arrange(comparison)

render_table <- function (data, .title) {
  last_row <- nrow(data)
  gt(data, rowname_col = "comparison") |>
    tab_header(title = md(.title)) |>
    fmt_number(columns = 2:7, decimals = 3) |>
    cols_align(align = "center") %>%
    tab_footnote(
      "ROC AUC: Receiver Operating Characteristic Area Under the Curve",
      locations = cells_column_labels(columns = `ROC AUC`)
    ) |>
    tab_footnote(
      "PR AUC: Precision-Recall Area Under the Curve",
      locations = cells_column_labels(columns = `PR AUC`)
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = list(cells_body(rows = last_row), cells_stub(rows = last_row))
    )
}

complex_table <- metrics %>%
  filter(model == "Fusion") %>%
  select(-model) %>%
  render_table("**HCC Fusion** Classifier")
simple_table <- metrics %>%
  filter(model == "Slim") %>%
  select(-model) %>%
  render_table("**HCC Slim** Classifier")
gtsave(complex_table, snakemake@output[[1]])
gtsave(simple_table, snakemake@output[[2]])
