library(tidyverse)
library(rstatix)

# cor_result <- read_csv("results/data/cor_with_clinical/glycan_cor_with_clinical.csv")
cor_result <- read_csv(snakemake@input[[1]])

plot_data <- cor_result %>%
  add_significance("p.adj") %>%
  mutate(p.adj.signif = if_else(p.adj.signif == "ns", NA, p.adj.signif))

p <- ggplot(plot_data, aes(clinical, glycan)) +
  geom_tile(aes(fill = cor)) +
  geom_text(aes(label = p.adj.signif), nudge_y = -0.2, size = 3) +
  scale_fill_gradient2(high = "#D26F32", mid = "white", low = "#275D87") +
  labs(x = "Clinical Variables", y = "Glycans", fill = "rho") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
# tgutil::ggpreview(width = 6, height = 10)
ggsave(snakemake@output[[1]], plot = p, width = 6, height = 10)