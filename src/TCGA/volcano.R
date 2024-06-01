library(tidyverse)
library(ggrepel)

# dea_result <- read_csv("results/data/TCGA/dea_results.csv")
# glyco_genes <- read_csv("data/glycogenes.csv") %>%
#   pull(gene_name)
dea_result <- read_csv(snakemake@input[[1]])
glyco_genes <- read_csv(snakemake@input[[2]]) %>%
  pull(gene_name)

plot_data <- dea_result %>%
  select(logFC, FDR, gene_name) %>%
  mutate(
    label = if_else(
      (gene_name %in% glyco_genes) & (abs(logFC) > 1) & (FDR < 0.01),
      gene_name, NA
    ),
    regulate = case_when(
      (logFC > log2(2)) & (FDR < 0.01) ~ "up",
      (logFC < -log2(2)) & (FDR < 0.01) ~ "down",
      .default = "no"
    ),
  )

ggplot(plot_data, aes(logFC, -log10(FDR))) +
  geom_point(aes(color = regulate), shape = 16, size = 0.5) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_label_repel(aes(label = label), max.overlaps = 15, nudge_y = 5) +
  labs(
    x = expression(paste(log[2], "FC")),
    y = expression(paste(-log[10], "FDR")),
  ) +
  guides(color = "none") +
  scale_color_manual(values = c(up = "#CD0000", down = "#27408B", no = "grey")) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )
# tgutil::ggpreview(width = 4, height = 4)
ggsave(snakemake@output[[1]], width = 4, height = 4)