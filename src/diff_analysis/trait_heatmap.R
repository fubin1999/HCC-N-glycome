library(tidyverse)
library(ComplexHeatmap)
library(circlize)

trait_data <- read_csv("results/data/prepared/filtered_derived_traits.csv")
groups <- read_csv("results/data/prepared/groups.csv")
post_hoc_result <- read_csv("results/data/diff_analysis/trait_post_hoc.csv")

traits_to_show <- post_hoc_result %>%
  distinct(trait) %>%
  pull(trait)

data <- trait_data %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") %>%
  filter(trait %in% traits_to_show) %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

plot_data <- data %>%
  group_by(trait) %>%
  mutate(value = as.double(scale(value))) %>%
  group_by(trait, group) %>%
  summarise(mean_value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_value) %>%
  mutate(trait_type = case_when(
    str_ends(trait, "F") ~ "Fucosylation",
    str_ends(trait, "S") ~ "Sialylation",
    str_ends(trait, "G") ~ "Galactosylation",
    str_ends(trait, "B") ~ "Bisection",
    .default = "Complexity"
  )) %>%
  mutate(trait_type = factor(trait_type, levels = c(
    "Complexity", "Bisection", "Galactosylation", "Fucosylation", "Sialylation"
  )))

mat <- plot_data %>%
  select(-trait_type) %>%
  column_to_rownames("trait") %>%
  as.matrix() %>%
  t()

col_split <- plot_data %>%
  pull(trait_type)

col_fun <- colorRamp2(c(-1, 0, 1), c("#275D87", "white", "#D26F32"))

ht <- Heatmap(
  mat,
  name = "Z-score",
  col = col_fun,
  width = unit(ncol(mat) / 2.5, "cm"),
  height = unit(nrow(mat) / 2.5, "cm"),
  cluster_rows = FALSE,
  column_names_rot = 60,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10),
  column_split = col_split,
  cluster_column_slices = FALSE,
  rect_gp = gpar(col = "white", lwd = 1.5)
)
ht <- draw(ht)
w <- ComplexHeatmap:::width(ht)
w <- convertX(w, "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(ht)
h <- convertY(h, "inch", valueOnly = TRUE)

pdf(snakemake@output[[1]], width = w, height = h)
draw(ht)
dev.off()

write_csv(plot_data, "results/source_data/Supplementary_Figure_8.csv")
