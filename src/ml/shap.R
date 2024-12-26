library(tidyverse)
library(ggbeeswarm)

shap <- read_csv("results/data/ml/shap.csv")
test_data <- read_csv("results/data/ml/test_data.csv")

importance <- shap %>% 
  pivot_longer(-c(model, sample), names_to = "glycan", values_to = "shap") %>% 
  summarise(importance = mean(abs(shap)), .by = c(model, glycan)) %>% 
  slice_max(importance, n = 15, by = model)

ranks <- test_data %>% 
  select(-group) %>% 
  pivot_longer(-sample, names_to = "glycan", values_to = "abundance") %>% 
  mutate(rank = rank(abundance), .by = glycan) %>% 
  select(-abundance)

plot_data <- shap %>% 
  pivot_longer(-c(model, sample), names_to = "glycan", values_to = "shap") %>% 
  left_join(ranks, by = c("sample", "glycan")) %>% 
  right_join(importance, by = c("model", "glycan"))

plot_shap_beeswarm <- function(data) {
  ggplot(data, aes(shap, reorder(glycan, importance), color = rank)) +
    geom_quasirandom() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_gradient(low = "#3A3D8F", high = "#ED6A5A") +
    labs(x = "SHAP value", y = NULL, color = "Feature Value") +
    theme_classic() +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
    )
}

plot_df <- plot_data %>% 
  nest_by(model) %>% 
  mutate(plot = list(plot_shap_beeswarm(data)))

merged_shap_p <- reduce(plot_df$plot, `+`) + 
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")
tgutil::ggpreview(width = 12, height = 4)
ggsave("results/figures/ml/shap_beeswarm.pdf", merged_shap_p, width = 12, height = 4)
