library(tidyverse)
library(yardstick)

# predictions <- read_csv("results/data/ml/predictions.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
predictions <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
predictions <- predictions |> 
  inner_join(groups, by = "sample") |> 
  mutate(
    target = factor(target, levels = c(T, F)), 
    prediction = factor(prediction, levels = c(T, F))
  )

complex_global_predictions <- predictions |> filter(model == "HCC Fusion")
complex_HC_predictions <- predictions |> filter(group %in% c("HC", "HCC"), model == "HCC Fusion")
complex_MC_predictions <- predictions |> filter(group %in% c("CHB", "HCC"), model == "HCC Fusion")
complex_YC_predictions <- predictions |> filter(group %in% c("LC", "HCC"), model == "HCC Fusion")
simple_global_predictions <- predictions |> filter(model == "HCC Slim")
simple_HC_predictions <- predictions |> filter(group %in% c("HC", "HCC"), model == "HCC Slim")
simple_MC_predictions <- predictions |> filter(group %in% c("CHB", "HCC"), model == "HCC Slim")
simple_YC_predictions <- predictions |> filter(group %in% c("LC", "HCC"), model == "HCC Slim")

get_metrics <- function(data) {
  class_metrics <- metric_set(accuracy, sensitivity, specificity, f_meas)
  proba_metrics <- metric_set(roc_auc, pr_auc)
  df1 <- class_metrics(data, truth = target, estimate = prediction)
  df2 <- proba_metrics(data, probability, truth = target)
  bind_rows(df1, df2)
}

data_list <- list(
  `Fusion_Control/HCC` = complex_global_predictions,
  `Fusion_HC/HCC` = complex_HC_predictions,
  `Fusion_CHB/HCC` = complex_MC_predictions,
  `Fusion_LC/HCC` = complex_YC_predictions,
  `Slim_Control/HCC` = simple_global_predictions,
  `Slim_HC/HCC` = simple_HC_predictions,
  `Slim_CHB/HCC` = simple_MC_predictions,
  `Slim_LC/HCC` = simple_YC_predictions
)
result <- map(data_list, get_metrics) |> 
  bind_rows(.id = "data") |>
  rename(metric = .metric, score = .estimate) |> 
  select(-.estimator) |> 
  sjmisc::rec(
    metric, 
    rec = "f_meas=f1_score;else=copy", 
    suffix = ""
  ) %>%
  separate_wider_delim(data, delim = "_", names = c("model", "comparison"))

write_csv(result, snakemake@output[[1]])