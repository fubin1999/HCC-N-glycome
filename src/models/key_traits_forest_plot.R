library(tidyverse)
library(patchwork)

# Load data-----
param_data <- read_csv("results/data/models/trait_parameters.csv")

param_data <- param_data %>%
  filter(trait %in% c("CA2", "CA3", "CA4", "CG", "TB", "TC", "TM", "THy", "TF", "TS")) %>%
  mutate(trait = case_match(
    trait,
    "CA2" ~ "Biantennary",
    "CA3" ~ "Triantennary",
    "CA4" ~ "Tetraantennary",
    "CG" ~ "Galactosylation",
    "TB" ~ "Bisecting",
    "TC" ~ "Complex Type",
    "TM" ~ "High-Mannose Type",
    "THy" ~ "Hybrid Type",
    "TF" ~ "Fucosylation",
    "TS" ~ "Sialylation"
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
  nest_by(trait) %>%
  mutate(plot = list(plot_forest(data) + ggtitle(trait))) %>%
  arrange(match(trait, c("Biantennary", "Triantennary", "Tetraantennary", "Bisecting",
                         "Complex Type", "High-Mannose Type", "Hybrid Type",
                         "Bisecting", "Galactosylation", "Fucosylation", "Sialylation")))

final_p <- reduce(plot_df$plot, `+`) + 
  plot_layout(nrow = 3) +
  plot_annotation(tag_levels = "a")
# tgutil::ggpreview(final_p, width = 15, height = 6)
ggsave("results/figures/models/key_trait_forest_plot.pdf", final_p, width = 12, height = 9)

source_data_df <- param_data %>%
  nest_by(trait) %>%
  ungroup() %>%
  mutate(
    tag = letters[1:10],
    filepath = str_glue("results/source_data/Supplementary_Figure_11{tag}.csv")
  )

walk2(source_data_df$data, source_data_df$filepath, ~write_csv(.x, .y))
