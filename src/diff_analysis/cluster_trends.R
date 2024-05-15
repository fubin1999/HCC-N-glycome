source("renv/activate.R")

library(tidyverse)
library(patchwork)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clusters <- read_rds("results/data/diff_analysis/glycan_clusters.rds")
abundance <- read_csv(snakemake@input[["abundance"]])
groups <- read_csv(snakemake@input[["groups"]])
clusters <- read_rds(snakemake@input[["clusters"]])

data <- abundance |> 
  inner_join(groups, by = "sample") |> 
  filter(group != "QC") |>
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) |> 
  mutate(value = log(value)) |> 
  group_by(glycan) |> 
  mutate(value = (value - mean(value)) / sd(value)) |> 
  ungroup() |> 
  summarise(median = median(value), .by = c(glycan, group))

plot_trend <- function(idx) {
  title <- str_glue("Cluster {names(clusters)[[idx]]}")
  data |> 
    filter(glycan %in% clusters[[idx]]) |> 
    ggplot(aes(group, median, group = glycan)) +
    labs(x = "", y = "z-score median", title = title) +
    geom_point(color = "steelblue") +
    geom_line(color = "steelblue") +
    theme_classic()
}

plot_trend(1) + plot_trend(2) + plot_trend(3) + plot_trend(4)

tgutil::ggpreview(width = 6, height = 6)
ggsave(snakemake@output[[1]], width = 6, height = 6)
