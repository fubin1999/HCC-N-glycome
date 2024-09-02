library(tidyverse)
library(yardstick)

# preds <- read_csv("results/data/cor_with_liver_function/ALBI_model_preds.csv")
# scores <- read_csv("results/data/cor_with_liver_function/ALBI_model_scores.csv")

preds <- read_csv(snakemake@input[[1]])
scores <- read.csv(snakemake@input[[2]])

data_for_plot <- bind_rows(list(
  preds %>% mutate(fold = as.character(fold)),
  preds %>% mutate(fold = "all"))
) %>%
  mutate(truth = as.factor(truth)) %>%
  group_by(fold) %>%
  roc_curve(truth, prob.I)

micro_auc <- scores %>%
  filter(average == "micro", metric == "ROC AUC") %>%
  pull(score)

macro_auc <- scores %>%
  filter(average == "macro", metric == "ROC AUC") %>%
  pull(score)

roc_p <- data_for_plot %>%
  mutate(group = if_else(fold == "all", "Micro Average", "Cross Validation")) %>%
  ggplot(aes(1 - specificity, sensitivity, group = fold)) +
  geom_path(aes(alpha = group), color = "#275D87") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  annotate("text", x = 0.7, y = 0.2, label = paste("Micro AUC:", round(micro_auc, 3))) +
  annotate("text", x = 0.7, y = 0.1, label = paste("Macro AUC:", round(macro_auc, 3))) +
  scale_alpha_manual(values = c("Micro Average" = 1, "Cross Validation" = 0.2)) +
  coord_equal() +
  labs(alpha = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank())
# tgutil::ggpreview(roc_p, width = 5, height = 3)

ggsave(snakemake@output[[1]], roc_p, width = 5, height = 3)