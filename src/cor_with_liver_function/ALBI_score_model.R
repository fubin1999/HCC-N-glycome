library(tidyverse)
library(mlr3verse)
library(patchwork)

# glycan_data <- read_csv("results/data/prepared/processed_abundance.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

glycan_data <- read_csv(snakemake@input[[1]])
clinical <- read_csv(snakemake@input[[2]])
groups <- read_csv(snakemake@input[[3]])

data <- glycan_data %>%
  right_join(clinical %>% select(sample, ALBI_score), by = "sample") %>%
  left_join(groups, by = "sample") %>%
  filter(ALBI_score < 0)

tsk <- as_task_regr(data, target = "ALBI_score")
tsk$set_col_roles("sample", roles = "name")
tsk$set_col_roles("group", roles = "stratum")

rf <- lrn("regr.ranger")
dummy <- lrn("regr.featureless")
cv10 <- rsmp("cv", folds = 10)
set.seed(123)
rr_rf <- resample(tsk, rf, cv10)
rr_dummy <- resample(tsk, dummy, cv10)

metrics_rf <- rr_rf$aggregate(msrs(c("regr.rmse", "regr.mae", "regr.rsq"), average = "micro"))
metrics_dummy <- rr_dummy$aggregate(msrs(c("regr.rmse", "regr.mae", "regr.rsq"), average = "micro"))

metrics_df <- bind_rows(
  list(
    rf = tibble(metric = names(metrics_rf), score = metrics_rf),
    dummy = tibble(metric = names(metrics_dummy), score = metrics_dummy)
  ),
  .id = "model"
) %>%
  mutate(metric = case_match(
    metric,
    "regr.rmse" ~ "RMSE",
    "regr.mae" ~ "MAE",
    "regr.rsq" ~ "R2"
  ))

write_csv(metrics_df, snakemake@output[[1]])

rf_preds <- rr_rf$prediction() %>% as.data.table() %>% as.tibble()
rf_p <- ggplot(rf_preds, aes(truth, response)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "True ALBI score", y = "Predicted ALBI score") +
  annotate("text", x = -1, y = -2.6, label = str_c("RMSE: ", round(metrics_rf[["regr.rmse"]], 3)), size = 5) +
  annotate("text", x = -1, y = -2.7, label = str_c("MAE: ", round(metrics_rf[["regr.mae"]], 3)), size = 5) +
  annotate("text", x = -1, y = -2.8, label = str_c("R2: ", round(metrics_rf[["regr.rsq"]], 3)), size = 5) +
  theme_minimal()

ggsave(snakemake@output[[2]], rf_p, width = 6, height = 6)