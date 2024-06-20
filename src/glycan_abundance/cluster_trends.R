library(tidyverse)
library(patchwork)

# eigen_glycans <- read_csv("results/data/glycan_abundance/eigen_glycans.csv")
# groups <- read_csv("results/data/prepared/groups.csv") %>%
#   mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

eigen_glycans <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]]) %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

data <- eigen_glycans %>%
  left_join(groups, by = "sample")

plot_data <- data %>%
  summarise(
    mean = mean(eigen_glycan),
    sd = sd(eigen_glycan),
    ci_upper = Rmisc::CI(eigen_glycan)[["upper"]],
    ci_lower = Rmisc::CI(eigen_glycan)[["lower"]],
    .by = c(group, cluster)
  )

plot_trend <- function (data, .title) {
  data %>%
    ggplot(aes(group, mean, group = 1)) +
    geom_ribbon(aes(ymin = ci_upper, ymax = ci_lower), alpha = 0.1, fill = "steelblue") +
    geom_line(color = "steelblue") +
    labs(title = .title, y = "Eigenglycan (scaled)") +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_line(),
      axis.title.x = element_blank()
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.01, 0.05)))
}

plots <- plot_data %>%
  nest_by(cluster) %>%
  mutate(
    title = str_c("Cluster ", cluster),
    plot = list(plot_trend(data, title))
  ) %>%
  select(cluster, plot)

final_p <- reduce(plots$plot, `+`) + plot_layout(nrow = 1)
# tgutil::ggpreview(plot = final_p, width = 12, height = 3)
ggsave(snakemake@output[[1]], width = 12, height = 3)