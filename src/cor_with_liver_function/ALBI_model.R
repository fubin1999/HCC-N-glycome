library(tidyverse)
library(mlr3verse)
library(data.table)

glycan_data <- read_csv("results/data/prepared/processed_abundance.csv")
clinical <- read_csv("results/data/prepared/clinical.csv")
groups <- read_csv("results/data/prepared/groups.csv")

glycan_data <- read_csv(snakemake@input[[1]])
clinical <- read_csv(snakemake@input[[2]])
groups <- read_csv(snakemake@input[[3]])

data <- glycan_data %>%
  right_join(clinical %>% select(sample, ALBI_stage), by = "sample") %>%
  mutate(ALBI_stage = if_else(ALBI_stage == "I", "I", "II/III")) %>%
  left_join(groups, by = "sample")

tsk <- as_task_classif(data, target = "ALBI_stage", id = "liver_function")
tsk$set_col_roles("group", "stratum")
tsk$set_col_roles("sample", "name")
rf <- lrn("classif.ranger", predict_type = "prob")
cv10 <- rsmp("cv", folds = 10)

set.seed(123)
rr <- resample(tsk, rf, cv10)

metric_names <- c("classif.auc", "classif.acc", "classif.fbeta",
                  "classif.prauc", "classif.precision", "classif.recall",
                  "classif.sensitivity", "classif.specificity")
micro_scores <- rr$aggregate(msrs(metric_names, average = "micro"))
macro_scores <- rr$aggregate(msrs(metric_names, average = "macro"))
scores_df <- bind_rows(list(
  micro = tibble(metric = names(micro_scores), score = micro_scores),
  macro = tibble(metric = names(macro_scores), score = macro_scores)
  ), .id = "average") %>%
  mutate(metric = case_match(
    metric,
    "classif.auc" ~ "ROC AUC",
    "classif.acc" ~ "Accuracy",
    "classif.fbeta" ~ "F1 Score",
    "classif.prauc" ~ "PR AUC",
    "classif.precision" ~ "Precision",
    "classif.recall" ~ "Recall",
    "classif.sensitivity" ~ "Sensitivity",
    "classif.specificity" ~ "Specificity"
  ))

preds <- rr$predictions()
preds_df <- bind_rows(map(preds, ~as.data.table(.x)), .id = "fold") %>%
  left_join(tsk$row_names, by = c("row_ids" = "row_id")) %>%
  select(row_id = row_ids, row_name, everything())

write_csv(scores_df, snakemake@output[[1]])
write_csv(preds_df, snakemake@output[[2]])