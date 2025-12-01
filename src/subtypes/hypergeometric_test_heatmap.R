library(tidyverse)
library(rstatix)
library(patchwork)

fisher_result <- read_csv("results/data/subtypes/subtype_with_categoric_clinical_fisher_result.csv")
post_hoc_result <- read_csv("results/data/subtypes/subtype_with_categoric_clinical_post_hoc_result.csv")

plot_heatmap <- function (data, x.title) {
  ggplot(data, aes(clinical_subtype, subtype, fill = -log10(p_value))) +
    geom_tile(linewidth = 1, color = "white") +
    geom_text(aes(label = p_value.signif)) +
    scale_fill_gradient(low = "#FFE8E2", high = "#ee6c4d", limits = c(0, 10)) +
    labs(x = x.title, y = "Glycan Subtype", fill = "-log10(p)") +
    coord_equal() +
    theme_minimal() +
    theme(panel.grid = element_blank())
}

plot_data <- post_hoc_result %>%
  mutate(subtype = paste0("S", subtype)) %>%
  add_significance("p_value") %>%
  mutate(p_value.signif = if_else(p_value.signif == "ns", "", p_value.signif))

plot_df <- plot_data %>%
  nest_by(clinical_variable) %>%
  mutate(
    x.title = case_match(
      clinical_variable,
      "child_pugh" ~ "Child-Pugh Class",
      "ALBI_stage" ~ "ALBI Stage",
      "TNM_stage" ~ "TNM Stage"
    ),
    x.title = factor(x.title, levels = c("Child-Pugh Class", "ALBI Stage", "TNM Stage")),
    plot = list(plot_heatmap(data, x.title))
  ) %>%
  arrange(x.title) %>%
  select(clinical_variable, plot)
p <- reduce(plot_df$plot, `+`) + plot_layout(nrow = 1, guides = "collect")
# tgutil::ggpreview(p, width = 8, height = 2.5)

ggsave(snakemake@output[[1]], p, width = 8, height = 2.5)
write_csv(plot_data, "results/source_data/Supplementary_Figure_15c.csv")
