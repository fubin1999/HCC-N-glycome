library(VennDiagram)
library(tidyverse)

# all_data <- read_csv("results/data/prepared/raw_abundance_full.csv")
# conf_data <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

all_data <- read_csv(snakemake@input[[1]])
conf_data <- read_csv(snakemake@input[[2]])
groups <- read_csv(snakemake@input[[3]])

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
conf_venn_data <- prepare_venn_data(conf_data)

plot_venn <- function (data) {
  venn.diagram(
    data,
    filename = NULL,
    disable.logging = TRUE,
    col = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A"),
  )
}

pdf(snakemake@output[[1]], width = 4, height = 4)
p1 <- plot_venn(all_venn_data)
grid.draw(p1)
dev.off()

pdf(snakemake@output[[2]], width = 4, height = 4)
p3 <- plot_venn(conf_venn_data)
grid.draw(p3)
dev.off()