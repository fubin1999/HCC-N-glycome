library(tidyverse)
library(patchwork)


perm_imp <- read_csv("results/data/ml/permutation_feature_importance.csv")

plot_perm_imp <- function(data) {
  data %>% 
    mutate(
      left_label = if_else(importance > 0, glycan, ""),
      right_label = if_else(importance <= 0, glycan, ""),
      left_signif = if_else(p_value < 0.05 & importance <= 0, "*", ""),
      right_signif = if_else(p_value < 0.05 & importance > 0, "*", ""),
    ) %>% 
    ggplot(aes(importance, reorder(glycan, importance))) +
    geom_col(aes(fill = importance)) +
    geom_vline(xintercept = 0, color = "black") +
    geom_text(aes(label = left_label, x = -0.001), hjust = 1, size = 3) +
    geom_text(aes(label = right_label, x = 0.001), hjust = 0, size = 3) +
    geom_text(aes(label = right_signif), nudge_y = -0.25, nudge_x = 0.002) +
    geom_text(aes(label = left_signif), nudge_y = -0.25, nudge_x = -0.002) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colorbar(title = "Importance", position = "inside")) +
    scale_x_continuous(expand = expansion(mult = c(0.3, 0.05))) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1) +
    theme_classic() +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      legend.position.inside = c(0.75, 0.3),
    )
}

plot_df <- perm_imp %>% 
  mutate(model = case_match(
    model, 
    "HCC_HC" ~ "HCC vs HC",
    "HCC_CHB" ~ "HCC vs CHB",
    "HCC_LC" ~ "HCC vs LC",
    "global" ~ "HCC vs Rest",
  )) %>% 
  mutate(model = factor(model, levels = c("HCC vs HC", "HCC vs CHB", "HCC vs LC", "HCC vs Rest"))) %>% 
  nest_by(model) %>% 
  arrange(model) %>% 
  mutate(plot = list(
    plot_perm_imp(data) + 
      ggtitle(model) +
      theme(plot.title = element_text(hjust = 0.5))
  ))

final_p <- reduce(plot_df$plot, `+`) + plot_layout(nrow = 1)
tgutil::ggpreview(final_p, width = 12, height = 4)
ggsave("results/figures/ml/permumtation_importance.pdf", final_p, width = 12, height = 4)
