library(tidyverse)

# mrmr_result <- read_csv("results/data/ml/mrmr_result.csv")
mrmr_result <- read_csv(snakemake@input[[1]])

p <- mrmr_result %>%
  slice_head(n = 30) %>%
  ggplot(aes(n_features)) +
  geom_ribbon(aes(ymin = score_mean - score_std, ymax = score_mean + score_std),
              alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = score_mean), color = "steelblue") +
  geom_point(aes(y = score_mean), color = "steelblue") +
  labs(x = "No. of features", y = "ROC AUC") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90")
  )
# tgutil::ggpreview(width = 4, height = 4)
ggsave(snakemake@output[[1]], plot = p, width = 4, height = 4)