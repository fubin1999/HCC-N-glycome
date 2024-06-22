library(VennDiagram)
library(tidyverse)

# full_raw <- read_csv("results/data/prepared/raw_abundance_full.csv")
# subset <- read_csv("results/data/prepared/raw_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

full_raw <- read_csv(snakemake@input[[1]])
subset <- read_csv(snakemake@input[[2]])
groups <- read_csv(snakemake@input[[3]])

full_raw <- full_raw %>%
  semi_join(subset, by = "sample")

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

full_venn_data <- prepare_venn_data(full_raw)
subset_venn_data <- prepare_venn_data(subset)

plot_venn <- function (data) {
  venn.diagram(
    data,
    filename = NULL,
    disable.logging = TRUE,
    fill = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A"),
  )
}

pdf(snakemake@output[[1]], width = 4, height = 4)
p1 <- plot_venn(full_venn_data)
grid.draw(p1)
dev.off()

pdf(snakemake@output[[2]], width = 4, height = 4)
p2 <- plot_venn(subset_venn_data)
grid.draw(p2)
dev.off()