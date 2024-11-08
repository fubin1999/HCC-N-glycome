library(tidyverse)
library(patchwork)

# Load data-----
# param_data <- read_csv("results/data/models/glycan_parameters.csv")
param_data <- read_csv(snakemake@input[[1]])

# Plot-----
plot_forest <- function(data) {
  ggplot(data, aes(Std_Coefficient, glycan)) +
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

HCC_p <- param_data %>% filter(Parameter == "Group: HCC") %>% plot_forest() + ggtitle("Group: HCC")
LC_p <- param_data %>% filter(Parameter == "Group: LC") %>% plot_forest() + ggtitle("Group: LC")
CHB_p <- param_data %>% filter(Parameter == "Group: CHB") %>% plot_forest() + ggtitle("Group: CHB")
ALBI_p <- param_data %>% filter(Parameter == "ALBI Score") %>% plot_forest() + ggtitle("ALBI Score")

final_p <- ALBI_p + CHB_p + LC_p + HCC_p + plot_layout(nrow = 1, axes = "collect_y")
# tgutil::ggpreview(final_p, width = 8, height = 8)
ggsave(snakemake@output[[1]], final_p, width = 8, height = 8)
