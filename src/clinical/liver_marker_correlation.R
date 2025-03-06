library(tidyverse)
library(corrplot)


clinical_data <- read.csv("results/data/prepared/clinical.csv")
liver_markers <- c("AST", "ALT", "GGT", "ALB", "TBIL", "TP", "AAR", "ALBI_score")

liver_marker_corr_mat <- clinical_data %>%
  select(all_of(liver_markers)) %>%
  cor(method = "spearman")

pdf("results/figures/clinical/liver_marker_correlation.pdf")
corrplot(
  liver_marker_corr_mat,
  method = "color",
  tl.col = "black",
  addgrid.col = "white"
)
dev.off()
