library(tidyverse)
library(ggrepel)

# post_hoc <- read_csv("results/data/diff_analysis/glycan_post_hoc.csv")
post_hoc <- read_csv(snakemake@input[[1]])

p <- post_hoc %>%
  mutate(
    logp = -log10(p.adj),
    comparison = str_c(group1, "vs", group2)
  ) %>%
  select(glycan, logp, comparison) %>%
  filter(comparison %in% c("HCvsLC", "HCvsHCC")) %>%
  pivot_wider(names_from = comparison, values_from = logp) %>%
  mutate(
    color = case_when(
      HCvsLC > -log10(0.05) & HCvsHCC > -log10(0.05) ~ "Both",
      HCvsLC > -log10(0.05) | HCvsHCC > -log10(0.05) ~ "One",
      .default = "Neither"
    ),
    color = factor(color, levels = c("Both", "One", "Neither")),
    label = if_else(color != "Neither", glycan, ""),
  ) %>%
  ggplot(aes(HCvsLC, HCvsHCC)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_point(aes(color = color), size = 4, alpha = 0.5) +
  geom_text_repel(
    aes(label = label, color = color),
    force = 10, nudge_x = 0.8, nudge_y = 0.1, min.segment.length = 0,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("Neither" = "grey", "Both" = "#CD0000", "One" = "#FF9494"),
    labels = c(
      "Neither" = "Neither",
      "Both" = "Up in Both",
      "One" = "Only up in HCC/LC"
    )
  )+
  labs(x = "-log10p (LC vs HC)", y = "-log10p (HCC vs HC)", color = "") +
  theme_classic()
# tgutil::ggpreview(plot = p, width = 6.5, height = 5)

ggsave(snakemake@output[[1]], plot = p, width = 6.5, height = 6.5)