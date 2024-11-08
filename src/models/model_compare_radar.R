library(tidyverse)
library(ggradar)
library(patchwork)

# Load data-----
# compare_data <- read_csv("results/data/models/model_comparison.csv")
compare_data <- read_csv(snakemake@input[[1]])

# Prepare data-----
plot_data <- compare_data %>%
  mutate(name = factor(name, levels = c("ALBI", "group", "group + ALBI", "group * ALBI"))) %>%
  select(gcm, name, AIC_wt, AICc_wt, BIC_wt, RMSE, Sigma, gcm, R2, R2_adjusted) %>%
  rename(AIC = AIC_wt, AICc = AICc_wt, BIC = BIC_wt, `R2(adj)` = R2_adjusted) %>%
  mutate(across(c(Sigma, RMSE), ~ -.))

# Plot-----
plot_radar <- function(data) {
  data %>%
    mutate(across(where(is.numeric), scales::rescale)) %>%
    ggradar(
      background.circle.colour = "white",
      values.radar = character(0)
    ) +
    ggsci::scale_color_npg()
}

suppressWarnings(
  plot_df <- plot_data %>%
    nest_by(gcm) %>%
    mutate(plot = list(
      plot_radar(data) +
        ggtitle(gcm) +
        theme(plot.title = element_text(hjust = 0.5, size = 15))
    ))
)

suppressWarnings(final_p <- reduce(plot_df$plot, `+`) + plot_layout(guides = "collect"))
# suppressWarnings(tgutil::ggpreview(final_p, width = 14, height = 8))
ggsave(snakemake@output[[1]], final_p, width = 14, height = 8)