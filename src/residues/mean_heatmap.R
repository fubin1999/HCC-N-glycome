library(tidyverse)

# residue_data <- read_csv("results/data/residues/glycan_residues.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

residue_data <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

p <- residue_data %>%
  pivot_longer(-sample, names_to = "residue", values_to = "n") %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  group_by(residue) %>%
  mutate(z_score = as.double(scale(n))) %>%
  group_by(residue, group) %>%
  summarise(mean_n = mean(n), mean_z_score = mean(z_score), .groups = "drop") %>%
  ggplot(aes(residue, group)) +
  geom_tile(aes(fill = mean_z_score), color = "grey30") +
  scale_fill_continuous(low = "#FFD8C8", high = "#FFA07A") +
  geom_text(aes(label = scales::number(mean_n, accuracy = 0.01))) +
  labs(x = "", y = "", fill = "Z-score") +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )
# tgutil::ggpreview(width = 4, height = 2)
ggsave(snakemake@output[[1]], plot = p, width = 4, height = 2)