library(tidyverse)

# trait_data <- read_csv("results/data/derived_traits/derived_traits.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
trait_data <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

data <- trait_data %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  group_by(trait) %>%
  mutate(value = as.vector(scale(value)))

traits_to_plot <- c("THy", "CA1", "CA2", "CA3", "CA4", "CFc", "CFa", "CS", "CG", "CB")

p <- data %>%
  filter(trait %in% traits_to_plot) %>%
  ggplot(aes(x = group, y = value, color = group)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1, shape = 16) +
  facet_wrap(~trait, scales = "free_y", nrow = 2) +
  labs(x = "", y = "Z-score") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "none"
  ) +
  scale_color_manual(values = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A"))
# tgutil::ggpreview(width = 10, height = 5)
ggsave(snakemake@output[[1]], plot = p, width = 10, height = 5)