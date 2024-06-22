library(tidyverse)
library(corrplot)

# eigen_glycans <- read_csv("results/data/glycan_abundance/eigen_glycans.csv")
eigen_glycans <- read_csv(snakemake@input[[1]])

wider <- eigen_glycans %>%
  pivot_wider(names_from = "cluster", values_from = "eigen_glycan", names_prefix = "GCM") %>%
  column_to_rownames("sample")

cor_mat <- cor(wider)
p_mat <- cor.mtest(wider)$p

pdf(snakemake@output[[1]], width = 3, height = 3)
corrplot(
  cor_mat,
  p.mat = p_mat,
  col = COL2('RdBu', 10),
  method = "color",
  type = "lower",
  diag = FALSE,
  cl.cex = 0.8,
  addCoef.col ='black',
  tl.col = "black",
  insig = "blank",
  addgrid.col = "grey"
)
dev.off()