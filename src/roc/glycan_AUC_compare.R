library(tidyverse)
library(patchwork)

# glycan_aucs <- read_csv("results/data/roc/glycan_auc.csv")
glycan_aucs <- read_csv(snakemake@input[[1]])

plot_bar <- function(data) {
  ggplot(data, aes(x = auc, y = reorder(glycan, auc))) +
    geom_col(aes(fill = auc)) +
    geom_text(aes(x = 0.02, label = glycan), hjust = 0, size = 2.5, color = "white") +
    geom_text(aes(label = sprintf("%.2f", auc), color = auc), hjust = -0.1, size = 2.5) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 1)) +
    scale_fill_distiller(palette = "YlOrRd", limit = c(0.5, 1)) +
    scale_color_distiller(palette = "YlOrRd", limit = c(0.5, 1)) +
    labs(x = "ROC AUC", y = NULL) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_text(size = 12)
    )
}

plot_df <- glycan_aucs %>% 
  nest_by(comparison) %>% 
  mutate(
    comparison = str_replace_all(comparison, "_", " "),
    comparison = factor(comparison, levels = c(
      "HCC vs HC", "HCC vs CHB", "HCC vs LC", "HCC vs rest"
    ))
  ) %>% 
  arrange(comparison) %>% 
  mutate(p = list(plot_bar(data) + ggtitle(comparison)))

p <- reduce(plot_df$p, `+`) + plot_layout(nrow = 1, guides = "collect")
# tgutil::ggpreview(p, width = 8, height = 8)
ggsave(snakemake@output[[1]], p, width = 8, height = 8)
