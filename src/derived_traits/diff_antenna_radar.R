library(tidyverse)
library(ggradar)
library(patchwork)

# Read data-----
# trait_data <- read_csv("results/data/derived_traits/derived_traits.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
trait_data <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

data <- trait_data %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

# Prepare data-----
traits_to_plot <- c(
  "A1Fc", "A2Fc", "A3Fc", "A4Fc",
  "A1Fa", "A2Fa", "A3Fa", "A4Fa",
  "A1B", "A2B", "A3B", "A4B",
  "A1S", "A2S", "A3S", "A4S",
  "A1G", "A2G", "A3G", "A4G"
)
unique_groups <- unique(data$group)

plot_data <- data %>%
  group_by(trait) %>%
  mutate(z_score = as.vector(scale(value))) %>%
  ungroup() %>%
  summarise(z_score = mean(z_score), .by = c(trait, group))

plot_data <- expand_grid(trait = traits_to_plot, group = unique_groups) %>%
  left_join(plot_data, by = c("trait", "group")) %>%
  mutate(z_score = ifelse(is.na(z_score), 0, z_score))

plot_data <- plot_data %>%
  separate_wider_regex(trait, c("A", ant = "\\d", trait = ".*")) %>%
  mutate(ant = as.integer(ant))

# Plot-----
plot_radar <- function (data) {
  data %>%
    select(-ant) %>%
    pivot_wider(names_from = trait, values_from = z_score) %>%
    ggradar(
      values.radar = c(-0.8, 0.0, 0.8),
      grid.min = -0.8,
      grid.mid = 0.0,
      grid.max = 0.8,
      background.circle.colour = "white",
      legend.position = "right",
      group.colours = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A"),
      fill = TRUE,
      fill.alpha = 0,
      gridline.mid.colour = "grey",
      gridline.max.linetype = "solid",
      gridline.mid.linetype = "dashed",
      gridline.min.linetype = "solid",
    )
}

ant1_p <- plot_radar(plot_data %>% filter(ant == 1))
ant2_p <- plot_radar(plot_data %>% filter(ant == 2))
ant3_p <- plot_radar(plot_data %>% filter(ant == 3))
ant4_p <- plot_radar(plot_data %>% filter(ant == 4))

ant1_p + ant2_p + ant3_p + ant4_p + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
# tgutil::ggpreview(width = 10, height = 10)
ggsave(snakemake@output[[1]], width = 8, height = 8)