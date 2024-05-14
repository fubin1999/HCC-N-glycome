source("renv/activate.R")

library(tidyverse)
library(yardstick)
library(patchwork)

# predictions <- read_csv("results/data/ml/predictions.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
predictions <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
predictions <- predictions |> 
  inner_join(groups, by = "sample")

get_cm <- function(data) {
  data |> 
    mutate(target = as.factor(target), prediction = as.factor(prediction)) |>
    conf_mat(truth = target, estimate = prediction)
}

plot_cm <- function(cm) {
  autoplot(cm, type = "heatmap") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      legend.position = 0
    ) +
    coord_equal() +
    scale_fill_gradient(low = "grey95", high = "steelblue", limits = c(0, 70)) +
    scale_x_discrete(labels = c("Control", "HCC")) +
    scale_y_discrete(labels = c("HCC", "Control"))
}

global_cm <- get_cm(predictions)
global_p <- plot_cm(global_cm) +
  ggtitle("Test Set: H+M+Y/C")

HC_cm <- get_cm(predictions |> filter(group %in% c("H", "C")))
HC_p <- plot_cm(HC_cm) +
  ggtitle("Test Set: H/C")

MC_cm <- get_cm(predictions |> filter(group %in% c("M", "C")))
MC_p <- plot_cm(MC_cm) +
  ggtitle("Test Set: M/C")

YC_cm <- get_cm(predictions |> filter(group %in% c("Y", "C")))
YC_p <- plot_cm(YC_cm) +
  ggtitle("Test Set: Y/C")

global_p | HC_p | MC_p | YC_p
ggsave(snakemake@output[[1]], width = 12, height = 3)
