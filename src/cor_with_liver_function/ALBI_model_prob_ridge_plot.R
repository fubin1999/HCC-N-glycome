library(tidyverse)
library(ggridges)

# preds <- read_csv("results/data/cor_with_liver_function/ALBI_model_preds.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

preds <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

plot_data <- preds %>%
  left_join(groups, by = c("row_name" = "sample")) %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

p <- ggplot(plot_data, aes(`prob.II/III`, group, fill = truth, color = truth, point_color = truth)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
  geom_density_ridges(
    jittered_points = TRUE, scale = .95,
    rel_min_height = .01, point_shape = "|", point_size = 3,
    position = position_points_jitter(height = 0)
  ) +
  labs(x = "Model Probability of ALBI II/III", fill = "ALBI stage") +
  theme_ridges(center = TRUE) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  guides(
    color = "none",
    point_color = "none",
    fill = guide_legend(
      override.aes = list(
        fill = c("#427296", "#D26F32"),
        color = NA, point_color = NA)
    )
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0))) +
  scale_fill_manual(values = c("#42729650", "#D26F3250")) +
  scale_color_manual(values = c("#427296", "#D26F32")) +
  scale_discrete_manual("point_color", values = c("#427296", "#D26F32"), guide = "none")

# tgutil::ggpreview(p, width = 5.5, height = 4)

ggsave(snakemake@output[[1]], p, width = 5.5, height = 4)