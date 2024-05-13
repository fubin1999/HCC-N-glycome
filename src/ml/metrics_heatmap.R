source("renv/activate.R")

library(tidyverse)

metrics_file <- snakemake@input[[1]]
metrics <- read_csv(metrics_file) |> 
  mutate(
    comparison = factor(comparison, levels = rev(c("H+M+Y/C", "H/C", "M/C", "Y/C"))),
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
  )

ggplot(metrics, aes(metric, comparison, fill = score)) +
  geom_tile(color = "grey30", linewidth = 0.4) +
  geom_text(aes(label = scales::number(score, accuracy = 0.001))) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = 0
  ) +
  scale_fill_gradient(low = "#FFC4C2", high = "#D66460") +
  scale_x_discrete(position = "top")
# tgutil::ggpreview(width = 6, height = 3)
ggsave(snakemake@output[[1]], width = 5, height = 2.5)

