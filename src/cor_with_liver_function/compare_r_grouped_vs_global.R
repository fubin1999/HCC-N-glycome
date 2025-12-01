library(tidyverse)
library(patchwork)

corr_result_global <- read_csv("results/data/cor_with_liver_function/global_cor_result_with_liver_functions.csv")
corr_result_grouped <- read_csv("results/data/cor_with_liver_function/grouped_cor_result_with_liver_functions.csv")

corr_result_global <- read_csv(snakemake@input[[1]])
corr_result_grouped <- read_csv(snakemake@input[[2]])

ALBI_result <- bind_rows(
  corr_result_global %>% mutate(group = "global"),
  corr_result_grouped
) %>%
  filter(feature_type == "glycan", clinical_variable == "ALBI_score")

plot_data <- ALBI_result %>%
  select(glycan = feature, cor, group) %>%
  pivot_wider(names_from = group, values_from = cor)

max_cor <- max(abs(ALBI_result$cor))

plot_func <- function(data) {
  ggplot(data, aes(global, grouped)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_point(aes(size = abs(global + grouped), fill = global + grouped), shape = 21) +
    scale_fill_gradient2(low = "#275D87", mid = "white", high = "#D26F32") +
    scale_size_continuous(range = c(1, 5)) +
    labs(x = "Corr. Coef. (Global)") +
    lims(x = c(-max_cor, max_cor), y = c(-max_cor, max_cor)) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

plot_list <- map(
  c("HC", "CHB", "LC", "HCC"),
  ~ plot_data %>%
    rename(grouped = all_of(.x)) %>%
    plot_func() +
    ylab(glue::glue("Corr. Coef. ({.x})"))
)
p <- reduce(plot_list, `+`) + plot_layout(nrow = 1)
# tgutil::ggpreview(width = 12, height = 3)
ggsave(snakemake@output[[1]], p, width = 12, height = 3)

parallel_plot_data <- ALBI_result %>% 
  filter(group != "global") %>%
  left_join(
    ALBI_result %>% filter(group == "global") %>% select(feature, cor), 
    by = "feature", suffix = c("", "_global")
  ) %>% 
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

ggplot(parallel_plot_data, aes(group, cor, color = cor_global)) +
  geom_line(aes(group = feature)) + 
  geom_point() +
  labs(x = "", y = "Spearman's rho\nin each group", color = "Global\nSpearman's\nrho") +
  scale_color_gradient2(low = "#275D87", mid = "white", high = "#D26F32") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )
# tgutil::ggpreview(width = 4, height = 3)
# ggsave("results/figures/cor_with_liver_function/ALBI_score_parallel_plot.pdf", width = 4, height = 3)

write_csv(parallel_plot_data, "results/source_data/Figure_2e.csv")
write_csv(plot_data, "results/source_data/Figure_2f.csv")
