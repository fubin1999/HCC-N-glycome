library(tidyverse)
library(patchwork)

# grouped_cor_result <- read_csv("results/data/cor_with_liver_function/grouped_cor_result_with_liver_functions.csv")
grouped_cor_result <- read_csv(snakemake@input[[1]])

plot_boxplot <- function (data, .title) {
  ggplot(data, aes(group, abs(cor))) +
    geom_boxplot(color = "#D26F32", fill = "#D26F32", alpha = 0.2) +
    #geom_jitter(color = "#D26F32", width = 0.3, alpha = 0.3) +
    labs(y = "| Spearman's rho |", title = .title) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank()
    )
}

plot_df <- grouped_cor_result %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  nest_by(clinical_variable) %>%
  mutate(title = str_c("Correlation with ", clinical_variable)) %>%
  mutate(plot = list(plot_boxplot(data, title)))

p <- reduce(plot_df$plot, `+`) &
  ylim(c(0, 0.65))
# tgutil::ggpreview(p, width = 8, height = 6)

ggsave(snakemake@output[[1]], plot = p, width = 8, height = 6)