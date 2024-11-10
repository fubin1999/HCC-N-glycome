library(tidyverse)
library(patchwork)

# Load data-----
# param_data <- read_csv("results/data/models/glycan_parameters.csv")
param_data <- read_csv(snakemake@input[[1]])
var_name <- snakemake@params[["var_name"]]

# Plot-----
plot_forest <- function(data) {
  ggplot(data, aes(Std_Coefficient, .data[[var_name]])) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.4) +
    geom_point(aes(fill = Std_Coefficient), color = "black", size = 2.5, shape = 22) +
    scale_fill_gradient2(low = "#275D87", mid = "white", high = "#BE5800") +
    labs(x = "Std. Coefficient") +
    theme_classic() +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 10),
      panel.grid.major.y = element_line(color = "grey95", size = 0.5)
    )
}

plot_data <- param_data %>%
  filter(Parameter %in% c("ALBI Score", "Group: CHB", "Group: LC", "Group: HCC"))

max_x <- max(plot_data$CI_high)
min_x <- min(plot_data$CI_low)

plot_list <- plot_data %>%
  nest_by(Parameter) %>%
  mutate(plot = list(plot_forest(data) + ggtitle(Parameter))) %>%
  select(Parameter, plot) %>%
  deframe()

final_p <- reduce(plot_list, `+`) +
  plot_layout(nrow = 1, axes = "collect_y") &
  xlim(min_x, max_x)
# tgutil::ggpreview(final_p, width = 8, height = 8)
ggsave(snakemake@output[[1]], final_p, width = 8, height = 8)
