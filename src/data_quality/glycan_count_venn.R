library(VennDiagram)
library(tidyverse)

all_data <- read_csv("results/data/prepared/raw_abundance_full.csv")
struc_data <- read_csv("results/data/prepared/raw_abundance.csv")
conf_data <- read_csv("results/data/prepared/processed_abundance.csv")
groups <- read_csv("results/data/prepared/groups.csv")

prepare_venn_data <- function (data) {
  data %>%
    pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
    left_join(groups, by = "sample") %>%
    filter(!is.na(value)) %>%
    select(group, glycan) %>%
    filter(group != "QC") %>%
    group_by(group) %>%
    summarise(glycan = list(glycan)) %>%
    mutate(group = factor(group, levels = c("HC", "HCC", "CHB", "LC"))) %>%
    arrange(group) %>%
    deframe()
}

all_venn_data <- prepare_venn_data(all_data)
struc_venn_data <- prepare_venn_data(struc_data)
conf_venn_data <- prepare_venn_data(conf_data)

plot_venn <- function (data) {
  venn.diagram(
    data,
    filename = NULL,
    disable.logging = TRUE,
    col = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A"),
  )
}

output_files <- snakemake@output
venn_data <- list(all_venn_data, struc_venn_data, conf_venn_data)

# 循环生成 PDF 文件
for (i in seq_along(output_files)) {
  pdf(output_files[[i]], width = 4, height = 4)
  p <- plot_venn(venn_data[[i]])
  grid.draw(p)
  dev.off()
}

prepare_source_data <- function(venn_data) {
  enframe(venn_data, "group", "glycan") %>%
    unnest(glycan)
}

all_venn_data %>%
  prepare_source_data() %>%
  write_csv("results/source_data/Supplementary_Figure_4b.csv")
struc_venn_data %>%
  prepare_source_data() %>%
  write_csv("results/source_data/Supplementary_Figure_4c.csv")
conf_venn_data %>%
  prepare_source_data() %>%
  write_csv("results/source_data/Supplementary_Figure_4d.csv")
