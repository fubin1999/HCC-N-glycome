library(tidyverse)
library(ggsignif)

# residues <- read_csv("results/data/residues/glycan_residues.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
residues <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

data <- residues %>%
  pivot_longer(-sample, names_to = "residue", values_to = "value") %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  sjmisc::rec(residue, rec = "F=Fucose;G=Galactose;M=Mannose;N=GlcNAc;S=NeuAc", suffix = "")

p <- data %>%
  ggplot(aes(group, value)) +
  geom_boxplot(aes(color = group)) +
  geom_jitter(aes(color = group), width = 0.3) +
  geom_signif(
    comparisons = list(c("HC", "CHB"), c("HC", "LC"), c("HC", "HCC")),
    step_increase = 0.12,
    map_signif_level = TRUE
  ) +
  facet_wrap(~ residue, scales = "free") +
  labs(y = "Scaled Number per Glycan") +
  scale_color_manual(values = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(color = "none")
# tgutil::ggpreview(width = 9, height = 6)
ggsave(snakemake@output[[1]], plot = p, width = 9, height = 6)