library(tidyverse)
library(patchwork)

# mp_table <- read_csv("results/data/prepared/meta_properties.csv")
mp_table <- read_csv(snakemake@input[[1]])

fill_color <- "steelblue"
alpha <- 0.1
nudge_y <- 2.5

sia_p <- mp_table %>%
  summarise(n = n(), .by = nS) %>%
  ggplot(aes(x = nS, y = n)) +
  geom_col(fill = fill_color, color = fill_color, alpha = alpha) +
  geom_text(aes(label = n), nudge_y = nudge_y) +
  labs(x = "No. of Sialic Acids", y = "No. of Glycans")
sia_width <- length(unique(mp_table$nS))
# tgutil::ggpreview(plot = sia_p, width = 3, height = 3)

fuc_p <- mp_table %>%
  summarise(n = n(), .by = nF) %>%
  ggplot(aes(x = nF, y = n)) +
  geom_col(fill = fill_color, color = fill_color, alpha = alpha) +
  geom_text(aes(label = n), nudge_y = nudge_y) +
  labs(x = "No. of Fucoses", y = "No. of Glycans")
fuc_width <- length(unique(mp_table$nF))

ant_p <- mp_table %>%
  summarise(n = n(), .by = nAnt) %>%
  filter(nAnt > 0) %>%
  ggplot(aes(x = nAnt, y = n)) +
  geom_col(fill = fill_color, color = fill_color, alpha = alpha) +
  geom_text(aes(label = n), nudge_y = nudge_y) +
  labs(x = "No. of Antennas", y = "No. of Glycans")
ant_width <- 4

type_p <- mp_table %>%
  summarise(n = n(), .by = type) %>%
  mutate(type = case_match(type, "complex" ~ "C", "high_mannose" ~ "M", "hybrid" ~ "H")) %>%
  ggplot(aes(x = type, y = n)) +
  geom_col(fill = fill_color, color = fill_color, alpha = alpha) +
  geom_text(aes(label = n), nudge_y = nudge_y) +
  labs(x = "Type of Glycan", y = "No. of Glycans")
type_width <- 3

bisect_p <- mp_table %>%
  summarise(n = n(), .by = B) %>%
  mutate(B = if_else(B, "Yes", "No")) %>%
  ggplot(aes(x = B, y = n)) +
  geom_col(fill = fill_color, color = fill_color, alpha = alpha) +
  geom_text(aes(label = n), nudge_y = nudge_y) +
  labs(x = "Bisecting", y = "No. of Glycans")
bisect_width <- 2

final_p <- (sia_p | fuc_p | ant_p | type_p | bisect_p) +
  plot_layout(
    widths = c(sia_width, fuc_width, ant_width, type_width, bisect_width),
    axes = "collect_y"
) &
  scale_y_continuous(
    limits = c(0, 55),
    expand = expansion(mult = c(0, 0.1))
  ) &
  guides(fill = "none") &
  theme_classic() &
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_text(hjust = 0)
  )
# tgutil::ggpreview(plot = final_p, width = 8, height = 3)

ggsave(snakemake@output[[1]], plot = final_p, width = 8, height = 3)