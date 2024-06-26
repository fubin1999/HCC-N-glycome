library(tidyverse)
library(patchwork)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# mp_table <- read_csv("results/data/prepared/meta_properties.csv")

abundance <- read_csv(snakemake@input[[1]])
mp_table <- read_csv(snakemake@input[[2]])

fill_color <- "steelblue"
alpha <- 0.1
nudge_y <- 2.5

sia_width <- length(unique(mp_table$nS))
fuc_width <- length(unique(mp_table$nF))
ant_width <- 4
type_width <- 3
bisect_width <- 2

sia_count_p <- mp_table %>%
  summarise(n = n(), .by = nS) %>%
  ggplot(aes(x = nS, y = n)) +
  geom_col(fill = fill_color, color = fill_color, alpha = alpha) +
  geom_text(aes(label = n), nudge_y = nudge_y) +
  labs(x = "No. of Sialic Acids", y = "No. of Glycans")

fuc_count_p <- mp_table %>%
  summarise(n = n(), .by = nF) %>%
  ggplot(aes(x = nF, y = n)) +
  geom_col(fill = fill_color, color = fill_color, alpha = alpha) +
  geom_text(aes(label = n), nudge_y = nudge_y) +
  labs(x = "No. of Fucoses", y = "No. of Glycans")

ant_count_p <- mp_table %>%
  summarise(n = n(), .by = nAnt) %>%
  filter(nAnt > 0) %>%
  ggplot(aes(x = nAnt, y = n)) +
  geom_col(fill = fill_color, color = fill_color, alpha = alpha) +
  geom_text(aes(label = n), nudge_y = nudge_y) +
  labs(x = "No. of Antennas", y = "No. of Glycans")

type_count_p <- mp_table %>%
  summarise(n = n(), .by = type) %>%
  mutate(type = case_match(type, "complex" ~ "C", "high_mannose" ~ "M", "hybrid" ~ "H")) %>%
  ggplot(aes(x = type, y = n)) +
  geom_col(fill = fill_color, color = fill_color, alpha = alpha) +
  geom_text(aes(label = n), nudge_y = nudge_y) +
  labs(x = "Type of Glycan", y = "No. of Glycans")

bisect_count_p <- mp_table %>%
  summarise(n = n(), .by = B) %>%
  mutate(B = if_else(B, "Yes", "No")) %>%
  ggplot(aes(x = B, y = n)) +
  geom_col(fill = fill_color, color = fill_color, alpha = alpha) +
  geom_text(aes(label = n), nudge_y = nudge_y) +
  labs(x = "Bisecting", y = "No. of Glycans")

final_count_p <- (sia_count_p | fuc_count_p | ant_count_p | type_count_p | bisect_count_p) +
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
# tgutil::ggpreview(plot = final_count_p, width = 8, height = 3)

abund_data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  group_by(sample) %>%
  mutate(value = value / sum(value)) %>%
  ungroup() %>%
  left_join(mp_table, by = "glycan")

sia_abund_p <- abund_data %>%
  summarise(value = sum(value), .by = c(nS, sample)) %>%
  mutate(nS = factor(nS)) %>%
  ggplot(aes(x = nS, y = value)) +
  geom_boxplot(color = fill_color, fill = fill_color, alpha = alpha) +
  xlab("No. of Sialic Acids")

fuc_abund_p <- abund_data %>%
  summarise(value = sum(value), .by = c(nF, sample)) %>%
  mutate(nF = factor(nF)) %>%
  ggplot(aes(x = nF, y = value)) +
  geom_boxplot(color = fill_color, fill = fill_color, alpha = alpha) +
  xlab("No. of Fucoses")

ant_abund_p <- abund_data %>%
  summarise(value = sum(value), .by = c(nAnt, sample)) %>%
  filter(nAnt > 0) %>%
  mutate(nAnt = factor(nAnt)) %>%
  ggplot(aes(x = nAnt, y = value)) +
  geom_boxplot(color = fill_color, fill = fill_color, alpha = alpha) +
  xlab("No. of Antennas")

type_abund_p <- abund_data %>%
  summarise(value = sum(value), .by = c(type, sample)) %>%
  mutate(type = case_match(type, "complex" ~ "C", "high_mannose" ~ "M", "hybrid" ~ "H")) %>%
  ggplot(aes(x = type, y = value)) +
  geom_boxplot(color = fill_color, fill = fill_color, alpha = alpha) +
  xlab("Type of Glycan")

bisect_abund_p <- abund_data %>%
  summarise(value = sum(value), .by = c(B, sample)) %>%
  mutate(B = if_else(B, "Yes", "No")) %>%
  ggplot(aes(x = B, y = value)) +
  geom_boxplot(color = fill_color, fill = fill_color, alpha = alpha) +
  xlab("Bisecting")

final_abund_p <- (sia_abund_p | fuc_abund_p | ant_abund_p | type_abund_p | bisect_abund_p) +
  plot_layout(
    widths = c(sia_width, fuc_width, ant_width, type_width, bisect_width),
    axes = "collect_y"
) &
  scale_y_continuous(
    name = "Relative Abundance",
    limits = c(0, 0.5),
    expand = expansion(mult = c(0.05, 0.05)),
    labels = scales::percent_format(accuracy = 1)
  ) &
  guides(fill = "none") &
  theme_classic()
# tgutil::ggpreview(plot = final_abund_p, width = 8, height = 3)

ggsave(snakemake@output[[1]], plot = final_count_p, width = 8, height = 3)
ggsave(snakemake@output[[2]], plot = final_abund_p, width = 8, height = 3)