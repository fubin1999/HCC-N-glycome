source("renv/activate.R")

library(tidyverse)
library(probably)

# predictions <- read_csv("results/data/ml/predictions.csv")
predictions <- read_csv(snakemake@input[[1]])

cal_plot_logistic(predictions, target, probability, conf_level = 0.9, include_rug = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(snakemake@output[[1]], width = 4, height = 4)