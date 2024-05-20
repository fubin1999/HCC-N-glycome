source("renv/activate.R")

library(tidyverse)
library(patchwork)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clusters <- read_csv("results/data/diff_analysis/glycan_clusters.csv")
abundance <- read_csv(snakemake@input[["abundance"]])
groups <- read_csv(snakemake@input[["groups"]])
clusters <- read_csv(snakemake@input[["clusters"]])

data <- abundance |> 
  inner_join(groups, by = "sample") |> 
  filter(group != "QC") |>
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) |> 
  mutate(value = log(value)) |> 
  group_by(glycan) |> 
  mutate(value = (value - mean(value)) / sd(value)) |> 
  ungroup() |> 
  summarise(median = median(value), .by = c(glycan, group)) |> 
  inner_join(clusters, by = "glycan")

plot_trend <- function(data, .cluster) {
  data |> 
    ggplot(aes(group, median, group = glycan)) +
    geom_point(color = "steelblue") +
    geom_line(color = "steelblue") +
    labs(x = "", y = "z-score median", title = str_glue("Cluster {.cluster}")) +
    theme_classic()
}

trend_plots <- data |> 
  nest(.by = cluster) |> 
  arrange(cluster) |> 
  mutate(plot = map2(data, cluster, plot_trend))
wrap_plots(trend_plots$plot, ncol = 1)

# tgutil::ggpreview(width = 6, height = 6)
ggsave(snakemake@output[[1]], width = 3, height = 12)
