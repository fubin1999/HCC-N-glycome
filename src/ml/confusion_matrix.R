library(tidyverse)
library(yardstick)
library(patchwork)

# predictions <- read_csv("results/data/ml/predictions.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
predictions <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
predictions <- predictions |> 
  inner_join(groups, by = "sample")

plot_cm <- function(data) {
  cm <- data |>
    mutate(target = as.factor(target), prediction = as.factor(prediction)) |>
    conf_mat(truth = target, estimate = prediction)
  autoplot(cm, type = "heatmap") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    coord_equal() +
    scale_fill_gradient(low = "grey95", high = "steelblue", limits = c(0, 70)) +
    scale_x_discrete(labels = c("Control", "HCC")) +
    scale_y_discrete(labels = c("HCC", "Control"))
}

complex_global_p <- predictions %>%
  filter(model == "HCC Fusion") %>%
  plot_cm() +
  ggtitle("Test Set: Control/HCC", subtitle = "HCC Fusion")
complex_HC_p <- predictions %>%
  filter(model == "HCC Fusion", group %in% c("HC", "HCC")) %>%
  plot_cm() +
  ggtitle("Test Set: HC/HCC", subtitle = "HCC Fusion")
complex_MC_p <- predictions %>%
  filter(model == "HCC Fusion", group %in% c("CHB", "HCC")) %>%
  plot_cm() +
  ggtitle("Test Set: CHB/HCC", subtitle = "HCC Fusion")
complex_YC_p <- predictions %>%
  filter(model == "HCC Fusion", group %in% c("LC", "HCC")) %>%
  plot_cm() +
  ggtitle("Test Set: LC/HCC", subtitle = "HCC Fusion")

simple_global_p <- predictions %>%
  filter(model == "HCC Slim") %>%
  plot_cm() +
  ggtitle("Test Set: Control/HCC", subtitle = "HCC Slim")
simple_HC_p <- predictions %>%
  filter(model == "HCC Slim", group %in% c("HC", "HCC")) %>%
  plot_cm() +
  ggtitle("Test Set: HC/HCC", subtitle = "HCC Slim")
simple_MC_p <- predictions %>%
  filter(model == "HCC Slim", group %in% c("CHB", "HCC")) %>%
  plot_cm() +
  ggtitle("Test Set: CHB/HCC", subtitle = "HCC Slim")
simple_YC_p <- predictions %>%
  filter(model == "HCC Slim", group %in% c("LC", "HCC")) %>%
  plot_cm() +
  ggtitle("Test Set: LC/HCC", subtitle = "HCC Slim")

p <- (complex_global_p | complex_HC_p | complex_MC_p | complex_YC_p) /
  (simple_global_p | simple_HC_p | simple_MC_p | simple_YC_p)
# tgutil::ggpreview(width = 12, height = 6)
ggsave(snakemake@output[[1]], plot = p, width = 12, height = 6)
