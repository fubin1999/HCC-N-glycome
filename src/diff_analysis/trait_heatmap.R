library(tidyverse)
library(ComplexHeatmap)
library(circlize)

trait_data <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
post_hoc_result <- read_csv(snakemake@input[[3]])

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
  as.matrix()

row_split <- plot_data %>%
  pull(trait_type)

col_fun <- colorRamp2(c(-1, 0, 1), c("#275D87", "white", "#D26F32"))
pdf(snakemake@output[[1]])
Heatmap(
  mat,
  name = "Z-score",
  col = col_fun,
  width = unit(ncol(mat) / 2.5, "cm"),
  height = unit(nrow(mat) / 2.5, "cm"),
  cluster_columns = FALSE,
  row_split = row_split,
  cluster_row_slices = FALSE,
  rect_gp = gpar(col = "white", lwd = 1)
)
dev.off()