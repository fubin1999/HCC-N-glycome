library(tidyverse)
library(patchwork)

# preds <- read_csv("results/data/cor_with_liver_function/ALBI_model_preds.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

preds <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

plot_data <- preds %>%
  left_join(groups, by = c("row_name" = "sample")) %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  mutate(
    truth = if_else(truth == "I", "ALBI I", "ALBI II/III"),
    truth = factor(truth, levels = c("ALBI I", "ALBI II/III")),
    response = if_else(response == "I", "ALBI I", "ALBI II/III"),
    response = factor(response, levels = c("ALBI I", "ALBI II/III"))
  ) %>%
  group_by(group, truth, response) %>%
  count() %>%
  ungroup() %>%
  complete(truth, response, group, fill = list(n = 0))

max_n <- max(plot_data$n)
min_n <- min(plot_data$n)

plot_confusion_matrix <- function (data) {
  ggplot(data, aes(truth, response)) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = n), size = 5) +
    coord_equal() +
    guides(fill = "none") +
    labs(x = "Truth", y = "Response") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
    ) +
    scale_fill_gradient(low = "#FFE5D4", high = "#D26F32", limits = c(min_n, max_n))
}

plot_df <- plot_data %>%
  nest_by(group) %>%
  mutate(plot = list(plot_confusion_matrix(data) + ggtitle(group)))
group_p <- reduce(plot_df$plot, `+`) +
  plot_layout(nrow = 2, axis_titles = "collect", axes = "collect") &
  theme(
    axis.text = element_blank(),
    axis.title = element_blank()
  )
# tgutil::ggpreview(group_p, width = 4, height = 4)

total_p <- plot_data %>%
  summarise(n = sum(n), .by = c("truth", "response")) %>%
  plot_confusion_matrix() +
  scale_fill_gradient(low = "#C8DCEC", high = "#427296") +
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5),
  )
# tgutil::ggpreview(total_p, width = 4, height = 4)

ggsave(snakemake@output[[1]], total_p, width = 3, height = 3)
ggsave(snakemake@output[[2]], group_p, width = 4, height = 4)