library(tidyverse)
library(ComplexHeatmap)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# raw_abundance <- read_csv("results/data/prepared/raw_abundance.csv")
# mp_table <- read_csv("results/data/prepared/meta_properties.csv")

abundance <- read_csv(snakemake@input[[1]])
raw_abundance <- read_csv(snakemake@input[[2]])
mp_table <- read_csv(snakemake@input[[3]])

mp_names <- setdiff(colnames(mp_table), "glycan")

mean_abundance <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  group_by(sample) %>%
  mutate(rela_abund = value / sum(value)) %>%
  group_by(glycan) %>%
  summarise(mean_rela_abund = mean(rela_abund))

detect_rates <- raw_abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  summarise(detect_rate = mean(!is.na(value)), .by = glycan)

data <- mp_table %>%
  mutate(
    `Complex Type` = type == "complex",
    `Bisecting` = B,
    `High Branching` = nAnt > 2,
    `Fucosylation` = nF > 0,
    `Sialylation` = nS > 0,
  ) %>%
  select(-all_of(mp_names)) %>%
  left_join(mean_abundance, by = "glycan") %>%
  left_join(detect_rates, by = "glycan") %>%
  arrange(desc(mean_rela_abund))

mat <- data %>%
  select(-mean_rela_abund, -detect_rate) %>%
  column_to_rownames("glycan") %>%
  t()
mat[mat == FALSE] <- "No"
mat[mat == TRUE] <- "Yes"

ha <- HeatmapAnnotation(
  `Rela. Abund. (%)` = anno_barplot(
    data$mean_rela_abund * 100,
    gp = gpar(fill = "#CAE2F5"),
    border = FALSE
  ),
  `Detect Rate (%)` = anno_barplot(
    data$detect_rate * 100,
    gp = gpar(fill = "#CAE2F5"),
    border = FALSE
  )
)

pdf(snakemake@output[[1]], width = 12, height = 3)
Heatmap(
  mat,
  name = "Feature Presence",
  col = c(Yes = "steelblue", No = "grey90"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_gap = unit(2, "mm"),
  rect_gp = gpar(col = "white", lwd = 2),
  top_annotation = ha,
)
dev.off()