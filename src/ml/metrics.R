source("renv/activate.R")

library(tidyverse)
library(yardstick)

predictions <- read_csv("results/data/ml/predictions.csv")
groups <- read_csv("results/data/prepared/groups.csv")
predictions <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
predictions <- predictions |> 
  inner_join(groups, by = "sample") |> 
  mutate(
    target = factor(target, levels = c(T, F)), 
    prediction = factor(prediction, levels = c(T, F))
  )

HC_predictions <- predictions |> filter(group %in% c("H", "C"))
MC_predictions <- predictions |> filter(group %in% c("M", "C"))
YC_predictions <- predictions |> filter(group %in% c("Y", "C"))

get_metrics <- function(data) {
  class_metrics <- metric_set(accuracy, sensitivity, specificity, f_meas)
  proba_metrics <- metric_set(roc_auc, pr_auc)
  df1 <- class_metrics(data, truth = target, estimate = prediction)
  df2 <- proba_metrics(data, probability, truth = target)
  bind_rows(df1, df2)
}

data_list <- list(
  `H+M+Y/C` = predictions,
  `H/C` = HC_predictions,
  `M/C` = MC_predictions,
  `Y/C` = YC_predictions
)
result <- map(data_list, get_metrics) |> 
  bind_rows(.id = "comparison") |> 
  rename(metric = .metric, score = .estimate) |> 
  select(-.estimator) |> 
  sjmisc::rec(
    metric, 
    rec = "f_meas=f1_score;else=copy", 
    suffix = ""
  )

write_csv(result, snakemake@output[[1]])