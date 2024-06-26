library(tidyverse)
library(ggrepel)

# fold_change <- read_csv("results/data/diff_analysis/glycan_fold_change.csv")

fold_change <- read_csv(snakemake@input[[1]])

p <- fold_change %>%
  mutate(
    logFC = log2(FC),
    comparison = str_c(group1, "vs", group2),
  ) %>%
  select(glycan, logFC, comparison) %>%
  filter(comparison %in% c("HCvsLC", "HCvsHCC")) %>%
  pivot_wider(names_from = comparison, values_from = logFC) %>%
  mutate(
    color = case_when(
      HCvsLC > log2(1.5) & HCvsHCC > log2(1.5) ~ "Both",
      HCvsLC > log2(1.5) ~ "HCvsLC",
      HCvsHCC > log2(1.5) ~ "HCvsHCC",
      .default = "Neither"
    ),
    both_label = if_else(color == "Both", glycan, ""),
    HCC_label = if_else(color == "HCvsHCC", glycan, ""),
    LC_label = if_else(color == "HCvsLC", glycan, "")
  ) %>%
  ggplot(aes(HCvsLC, HCvsHCC)) +
  geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = log2(1.5), linetype = "dashed", color = "grey") +
  geom_point(aes(color = color), size = 4, alpha = 0.5) +
  geom_text_repel(
    aes(label = both_label, color = color),
    force = 10, nudge_x = 0.8, nudge_y = 0.1, min.segment.length = 0,
    show.legend = FALSE
  ) +
  geom_text_repel(
    aes(label = HCC_label, color = color),
    force = 10, nudge_x = -0.5, nudge_y = 0.1, min.segment.length = 0,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("Neither" = "grey", "Both" = "#CD0000", "HCvsHCC" = "#FF9494", "HCvsLC" = "#FF9494"),
    labels = c(
      "Neither" = "Neither",
      "Both" = "Up in Both",
      "HCvsHCC" = "Only Up in HCC",
      "HCvsLC" = "Only Up in LC"
    )
  )+
  labs(x = "log2FC (LC vs HC)", y = "log2FC (HCC vs HC)", color = "") +
  theme_classic()
# tgutil::ggpreview(plot = p, width = 4.5, height = 3)
ggsave(snakemake@output[[1]], plot = p, width = 4.5, height = 3)