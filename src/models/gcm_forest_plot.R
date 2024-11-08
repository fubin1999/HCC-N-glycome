library(tidyverse)
library(patchwork)

# Load data-----
# param_data <- read_csv("results/data/models/parameters.csv")
param_data <- read_csv(snakemake@input[[1]])

param_data <- param_data %>%
  filter(formula_name == "sex + age + group + ALBI") %>%
  mutate(Parameter = case_match(
    Parameter,
    "sexF" ~ "Sex: F",
    "groupLC" ~ "Group: LC",
    "groupHCC" ~ "Group: HCC",
    "groupCHB" ~ "Group: CHB",
    "ALBI_score" ~ "ALBI Score",
    "age" ~ "Age",
    "(Intercept)" ~ "Intercept"
  )) %>%
  mutate(Parameter = factor(
    Parameter, levels = rev(c("Group: HCC", "Group: LC", "Group: CHB", "ALBI Score",
                              "Sex: F", "Age", "Intercept"))
  ))

# Plot-----
plot_forest <- function(data) {
  ggplot(data, aes(Std_Coefficient, Parameter)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.1) +
    geom_point(aes(fill = Std_Coefficient), color = "black", size = 3, shape = 22) +
    scale_fill_gradient2(low = "#275D87", mid = "white", high = "#BE5800") +
    labs(x = "Std. Coefficient") +
    theme_classic() +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
}

plot_df <- param_data %>%
  nest_by(gcm) %>%
  mutate(plot = list(plot_forest(data) + ggtitle(gcm)))
final_p <- reduce(plot_df$plot, `+`) + plot_layout(nrow = 1)
# tgutil::ggpreview(final_p, width = 15, height = 3)
ggsave(snakemake@output[[1]], final_p, width = 15, height = 3)