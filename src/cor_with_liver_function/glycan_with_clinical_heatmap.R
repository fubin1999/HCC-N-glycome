library(tidyverse)
library(patchwork)

# cor_result <- read_csv("results/data/cor_with_liver_function/global_cor_result_with_liver_functions.csv")
# ttest_result <- read_csv("results/data/cor_with_liver_function/global_ttest_result_with_liver_functions.csv")

var_name <- snakemake@params[["var_name"]]
cor_result <- read_csv(snakemake@input[[1]]) %>% filter(feature_type == var_name)
ttest_result <- read_csv(snakemake@input[[2]]) %>% filter(feature_type == var_name)

numeric_p <- cor_result %>%
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

categoric_p <- ttest_result %>%
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

p <- numeric_p / categoric_p
ggsave(snakemake@output[[1]], p, width = 10, height = 4)
