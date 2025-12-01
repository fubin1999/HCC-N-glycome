library(tidyverse)
library(patchwork)
library(cowplot)

cor_result <- read_csv("results/data/cor_with_liver_function/grouped_cor_result_with_liver_functions.csv")
ttest_result <- read_csv("results/data/cor_with_liver_function/grouped_ttest_result_with_liver_functions.csv")

plot_numeric_heatmap <- function(cor_result) {
  cor_result %>% 
    mutate(relation = case_when(
      p.adj >= 0.05 ~ "Not Sig.",
      cor > 0 ~ "Pos.",
      cor < 0 ~ "Neg."
    )) %>%
    mutate(relation = factor(relation, levels = c("Pos.", "Neg.", "Not Sig."))) %>%
    ggplot(aes(feature, clinical_variable)) +
    geom_point(aes(size = abs(cor), color = relation)) +
    coord_equal() +
    labs(y = "Numeric", size = "| Spearman's rho |", color = "Relation") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "top"
    ) +
    scale_size_continuous(range = c(1, 3.5)) +
    scale_color_manual(values = c("Pos." = "#D26F32", "Neg." = "#275D87", "Not Sig." = "grey"))
}

plot_categoric_heatmap <- function(ttest_result) {
  ttest_result %>%
    mutate(
      regulate = case_when(
        p.adj >= 0.05 ~ "Not Sig.",
        effsize > 0 ~ "Up",
        effsize < 0 ~ "Down"
      ),
      regulate = factor(regulate, levels = c("Up", "Down", "Not Sig.")),
      clinical_variable = case_match(
        clinical_variable,
        "child_pugh" ~ "Child Pugh",
        "ALBI_stage" ~ "ALBI stage"
      )
    ) %>%
    ggplot(aes(feature, clinical_variable)) +
    geom_point(aes(size = abs(effsize), color = regulate)) +
    coord_equal() +
    labs(y = "Categoric", size = "| Cohen's d |", color = "Regulate") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
      axis.title.x = element_blank(),
      legend.position = "top"
    ) +
    scale_size_continuous(range = c(1, 3.5)) +
    scale_color_manual(values = c("Up" = "#8318B0", "Down" = "#429440", "Not Sig." = "grey"))
}

plot_merged_heatmap <- function(cor_result, ttest_result) {
  numeric_p <- plot_numeric_heatmap(cor_result)
  categoric_p <- plot_categoric_heatmap(ttest_result)
  numeric_p / categoric_p
}

LC_cor_result <- cor_result %>% 
  filter(feature_type == "glycan", group == "LC")
HCC_cor_result <- cor_result %>% 
  filter(feature_type == "glycan", group == "HCC")
LC_ttest_result <- ttest_result %>% 
  filter(feature_type == "glycan", group == "LC")
HCC_ttest_result <- ttest_result %>%
  filter(feature_type == "glycan", group == "HCC")

LC_heatmap <- plot_merged_heatmap(LC_cor_result, LC_ttest_result)
HCC_heatmap <- plot_merged_heatmap(HCC_cor_result, HCC_ttest_result)

p <- plot_grid(LC_heatmap, HCC_heatmap, ncol = 1, labels = "auto")
tgutil::ggpreview(p, width = 10, height = 8)
ggsave("results/figures/cor_with_liver_function/grouped_lft_heatmap.pdf", p, width = 10, height = 8)

write_csv(LC_cor_result, "results/source_data/Supplementary_Figure_6a_1.csv")
write_csv(LC_ttest_result, "results/source_data/Supplementary_Figure_6a_2.csv")
write_csv(HCC_cor_result, "results/source_data/Supplementary_Figure_6b_1.csv")
write_csv(HCC_ttest_result, "results/source_data/Supplementary_Figure_6b_2.csv")
