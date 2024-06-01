library(tidyverse)
library(ComplexHeatmap)

# clinical <- read_csv("results/data/prepared/clinical.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
clinical <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

data <- groups %>%
  inner_join(clinical, by = "sample") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  arrange(group) %>%
  column_to_rownames("sample")

pdf(snakemake@output[[1]], width = 12, height = 8)
set.seed(2)
ha <- HeatmapAnnotation(
  df = data,
  col = list(
    group = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A")
  )
)
ht <- Heatmap(
  matrix(data = NA, nrow = 0, ncol = nrow(data)),
  top_annotation = ha,
)
draw(ht, annotation_legend_side = "bottom")
dev.off()