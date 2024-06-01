library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# post_hoc_result <- read_csv("results/data/diff_analysis/posthoc_for_glycans.csv")
post_hoc_result <- read_csv(snakemake@input[[1]])

plot_data <- post_hoc_result %>%
  mutate(comparison = str_glue("{group1} / {group2}")) %>%
  mutate(signif = as.integer(p.adj < 0.05)) %>%
  select(glycan, comparison, signif) %>%
  pivot_wider(names_from = comparison, values_from = signif) %>%
  column_to_rownames("glycan") %>%
  as.matrix()

m <- make_comb_mat(plot_data)
pdf(snakemake@output[[1]], width = 8, height = 3)
UpSet(
  m,
  top_annotation = upset_top_annotation(m, annotation_name_rot = 90),
  right_annotation = NULL,
  comb_col = "#738CA6"
)
dev.off()