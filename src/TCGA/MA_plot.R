library(tidyverse)
library(ggrepel)

# dea_result <- read_csv("results/data/TCGA/dea_results.csv")
# glycogenes <- read_csv("data/glycogenes.csv")$gene_name

dea_result <- read_csv(snakemake@input[[1]])
glycogenes <- read_csv(snakemake@input[[2]])$gene_name

p <- dea_result %>%
  select(gene_name, logFC, logCPM, FDR) %>%
  mutate(
    label = if_else(
      (gene_name %in% glycogenes) & (abs(logFC) > 1),
      gene_name, NA
    ),
    regulate = case_when(
      (logFC > 1) & (FDR < 0.01) ~ "up",
      (logFC < -1) & (FDR < 0.01) ~ "down",
      .default = "no"
    )
  ) %>%
  ggplot(aes(logCPM, logFC)) +
  geom_point(aes(color = regulate), shape = 16, size = 0.5) +
  geom_hline(yintercept = c(1, -1), linetype = "dashed") +
  geom_label_repel(aes(label = label), nudge_x = 2, max.overlaps = 15) +
  labs(
    x = expression(paste(log[2], "CPM")),
    y = expression(paste(log[2], "FC")),
    color = "Regulate",
  ) +
  guides(color = "none") +
  scale_color_manual(values = c(up = "#CD0000", down = "#27408B", no = "grey")) +
  theme_classic()
ggsave(snakemake@output[[1]], plot = p, width = 4, height = 4)