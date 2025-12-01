library(tidyverse)
library(patchwork)

eigen_glycan <- read_csv("results/data/glycan_coexpr/eigen_glycans.csv")
subtypes <- read_csv("results/data/subtypes/consensus_cluster_result.csv")

plot_data <- eigen_glycan %>%
  mutate(gcm = paste0("GCM", cluster), .keep = "unused") %>%
  right_join(subtypes %>% rename(subtype = class), by = "sample") %>%
  mutate(subtype = paste0("S", subtype)) %>%
  summarise(mean = mean(eigen_glycan), .by = c(subtype, gcm))

subtype_p <- plot_data %>%
  ggplot(aes(subtype, reorder(gcm, desc(gcm)), fill = mean)) +
  geom_tile(linewidth = 1, color = "white") +
  geom_text(aes(label = scales::number(mean, accuracy = 0.01))) +
  scale_fill_gradient2(low = "#3d5a80", mid = "white", high = "#ee6c4d", midpoint = 0) +
  labs(
    x = "Molecular Subtype",
    y = "Glycan Coexpression Module",
    fill = "Mean Eigenvalue",
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
  )
# tgutil::ggpreview(subtype_p, width = 4, height = 3)

ggsave(snakemake@output[[1]], subtype_p, width = 4, height = 3)

write_csv(plot_data, "results/source_data/Figure_4c.csv")
