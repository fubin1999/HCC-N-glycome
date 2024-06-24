library(tidyverse)
library(ComplexHeatmap)
library(circlize)


# traits <- read_csv("results/data/prepared/filtered_derived_traits.csv")
# subtypes <- read_csv("results/data/subtype/cc_result.csv")
# diff_result <- read_csv("results/data/subtype/subtype_trait_diff.csv")

traits <- read_csv(snakemake@input[[1]])
subtypes <- read_csv(snakemake@input[[2]])
diff_result <- read_csv(snakemake@input[[3]])


diff_traits <- diff_result %>%
  filter(p.adj < 0.05) %>%
  pull(trait)

data <- traits %>%
  select(sample, all_of(diff_traits)) %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") %>%
  right_join(subtypes, by = "sample") %>%
  mutate(
    class = str_c("HCC_S", class),
    class = factor(class, levels = c("HCC_S1", "HCC_S2"))
  )

z_scores <- data %>%
  group_by(trait) %>%
  mutate(z_score = as.double(scale(value))) %>%
  group_by(trait, class) %>%
  summarise(mean_z_score = mean(z_score), .groups = "drop")

mat <- z_scores %>%
  pivot_wider(names_from = class, values_from = mean_z_score) %>%
  column_to_rownames("trait") %>%
  as.matrix()

col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("#4DBBD5FF", "white", "#E64B35FF"))
pdf(snakemake@output[[1]], width = 3, height = 6)
Heatmap(
  mat,
  name = "z-score",
  col = col_fun,
  cluster_columns = FALSE,
  rect_gp = gpar(col = "white", lwd = 2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10, col = "white"))
  }
)
dev.off()