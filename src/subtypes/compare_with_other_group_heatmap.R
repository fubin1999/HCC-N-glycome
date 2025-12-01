library(tidyverse)
library(ComplexHeatmap)
library(circlize)

abundance <- read_csv("results/data/prepared/processed_abundance.csv")
groups <- read_csv("results/data/prepared/groups.csv")
subtypes <- read_csv("results/data/subtypes/consensus_cluster_result.csv")

plot_data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  left_join(subtypes %>% mutate(subtype = paste0("S", class), .keep = "unused"), by = "sample") %>%
  mutate(group = if_else(group == "HCC", subtype, group)) %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "S1", "S2", "S3"))) %>%
  group_by(glycan) %>%
  mutate(z_score = log(value) %>% scale() %>% as.double()) %>%
  group_by(glycan, group) %>%
  summarize(mean_z_score = mean(z_score), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_z_score)

mat <- plot_data %>%
  column_to_rownames("glycan") %>%
  as.matrix()

col <- colorRamp2(c(-1.5, 0, 1.5), c("#3d5a80", "white", "#ee6c4d"))

ht <- Heatmap(
  mat,
  name = "z-score",
  col = col,
  cluster_columns = FALSE,
  rect_gp = gpar(col = "white", lwd = 2),
  column_split = c(rep("Ctrl.", 3), rep("HCC", 3)),
  column_gap = unit(2, "mm"),
  width = ncol(mat) * unit(4, "mm") + unit(2, "mm"),
  height = nrow(mat) * unit(4, "mm")
)
ht <- draw(ht)
w <- ComplexHeatmap:::width(ht)
w <- convertX(w, "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(ht)
h <- convertY(h, "inch", valueOnly = TRUE)

pdf(snakemake@output[[1]], width = w, height = h)
draw(ht)
dev.off()

write_csv(plot_data, "results/source_data/Supplementary_Figure_14.csv")
