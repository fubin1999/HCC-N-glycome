library(tidyverse)
library(ggsignif)
library(patchwork)

# pred <- read_csv("results/data/cor_with_clinical/liver_function_model_pred.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

pred <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

data <- pred %>%
  left_join(groups, by = "sample") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

plot_boxplots <- function (data, .target) {
  ggplot(data, aes(group, pred)) +
    geom_boxplot(aes(color = group)) +
    geom_jitter(aes(color = group), width = 0.25, alpha = 0.5) +
    scale_color_manual(values = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A")) +
    geom_signif(
      comparisons = list(c("HC", "CHB"), c("HC", "LC"), c("HC", "HCC")),
      map_signif_level = TRUE,
      step_increase = 0.1, vjust = 0
    ) +
    guides(color = FALSE) +
    labs(y = "Predicted value", title = str_glue("Predicted {.target}")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
}

plot_df <- data %>%
  nest_by(target) %>%
  mutate(plot = list(plot_boxplots(data, target))) %>%
  select(target, plot)

p <- reduce(plot_df$plot, `+`) +
  plot_layout(nrow = 1)
# tgutil::ggpreview(width = 12, height = 4)
ggsave(snakemake@output[[1]], p, width = 10, height = 3.5)
