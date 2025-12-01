library(tidyverse)
library(patchwork)
library(rstatix)
library(ggprism)

eigen_glycans <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]]) %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))
post_hoc_result <- read_csv(snakemake@input[[3]])

eigen_glycans <- read_csv("results/data/glycan_coexpr/eigen_glycans.csv")
groups <- read_csv("results/data/prepared/groups.csv") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))
post_hoc_result <- read_csv("results/data/glycan_coexpr/cluster_post_hoc.csv")

data <- eigen_glycans %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC")

plot_data <- data %>%
  summarise(
    mean = mean(eigen_glycan),
    sd = sd(eigen_glycan),
    ci_upper = Rmisc::CI(eigen_glycan)[["upper"]],
    ci_lower = Rmisc::CI(eigen_glycan)[["lower"]],
    .by = c(group, cluster)
  )

max_y_pos <- plot_data %>%
  group_by(cluster) %>%
  summarise(
    range = max(ci_upper) - min(ci_lower),
    max = max(ci_upper)
  )

p_value_data <- post_hoc_result %>%
  mutate(label = scales::scientific(p.adj)) %>%
  select(cluster, group1, group2, label) %>%
  left_join(max_y_pos, by = "cluster") %>%
  mutate(y.position = case_when(
    (group1 == "HC") & (group2 == "CHB") ~ max + range * 0.1,
    (group1 == "CHB") & (group2 == "LC") ~ max + range * 0.1,
    (group1 == "LC") & (group2 == "HCC") ~ max + range * 0.1,
    (group1 == "HC") & (group2 == "LC") ~ max + range * 0.2,
    (group1 == "CHB") & (group2 == "HCC") ~ max + range * 0.35,
    (group1 == "HC") & (group2 == "HCC") ~ max + range * 0.45
  )) %>%
  select(cluster, group1, group2, label, y.position)

plot_trend <- function (data, p_val_df, .title) {
  data %>%
    ggplot(aes(group, mean, group = 1)) +
    geom_ribbon(aes(ymin = ci_upper, ymax = ci_lower), alpha = 0.1, fill = "steelblue") +
    geom_line(color = "steelblue") +
    add_pvalue(p_val_df, bracket.size = 0.3) +
    labs(title = .title, y = "Eigenglycan (scaled)") +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_line(color = "gray90"),
      axis.title.x = element_blank()
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.01, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
}

plots <- plot_data %>%
  nest_by(cluster, .key = "value_data") %>%
  left_join(nest_by(p_value_data, cluster, .key = "p_data"), by = "cluster") %>%
  mutate(
    title = str_c("GCM", cluster),
    plot = list(plot_trend(value_data, p_data, title))
  ) %>%
  select(cluster, plot)

final_p <- reduce(plots$plot, `+`) + plot_layout(nrow = 1)
# tgutil::ggpreview(plot = final_p, width = 12, height = 3)
ggsave("results/figures/glycan_coexpr/glycan_cluster_trends.pdf", width = 12, height = 3)

write_csv(plot_data, "results/source_data/Figure_3f_1.csv")
write_csv(p_value_data, "results/source_data/Figure_3f_2.csv")
