library(tidyverse)
library(ggprism)
library(patchwork)

trait_data <- read_csv("results/data/prepared/filtered_derived_traits.csv")
groups <- read_csv("results/data/prepared/groups.csv")
post_hoc_result <- read_csv("results/data/diff_analysis/trait_post_hoc.csv")

col <- c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A")

data <- trait_data %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

trait_df <- tribble(
  ~trait, ~description, ~percent,
  "CA2", "Proportion of bi-antennary glycans\nwithin complex glycans", TRUE,
  "CA3", "Proportion of tri-antennary glycans\nwithin complex glycans", TRUE,
  "CA4", "Proportion of tetra-antennary glycans\nwithin complex glycans", TRUE,
  "CGS", "Average number of sialic acids per galactose\non complex glycans", FALSE,
  "CG", "Average number of galactoses per antenna\non complex glycans", FALSE,
  "TF", "Proportion of core-fucosylated glycans", TRUE,
  "TB", "Proportion of bisecting glycans", TRUE
)

plot_data <- data %>%
  nest_by(trait) %>%
  right_join(trait_df, by = "trait")

plot_boxplot <- function (data, .trait, .descr, .perc) {
  sub_post_hoc_result <- post_hoc_result %>%
    filter(trait == .trait)
  HC_CHB_p_val <- sub_post_hoc_result %>%
    filter(group1 == "HC", group2 == "CHB") %>%
    pull(p.adj)
  HC_LC_p_val <- sub_post_hoc_result %>%
    filter(group1 == "HC", group2 == "LC") %>%
    pull(p.adj)
  HC_HCC_p_val <- sub_post_hoc_result %>%
    filter(group1 == "HC", group2 == "HCC") %>%
    pull(p.adj)
  .max <- max(data$value)
  .min <- min(data$value)
  .range <- .max - .min
  p_val_df <- tribble(
    ~group1, ~group2, ~label, ~y.position,
    "HC", "CHB", HC_CHB_p_val, .max + 0.1 * .range,
    "HC", "LC", HC_LC_p_val, .max + 0.25 * .range,
    "HC", "HCC", HC_HCC_p_val, .max + 0.4 * .range,
  ) |>
    mutate(label = scales::scientific(label))
  p <- ggplot(data, aes(group, value)) +
    geom_boxplot(aes(color = group), outlier.alpha = 0) +
    geom_jitter(aes(color = group), alpha = 0.5, width = 0.3, shape = 16, size = 1) +
    scale_color_manual(values = col) +
    add_pvalue(p_val_df, bracket.size = 0.3) +
    guides(color = "none") +
    labs(
      title = .trait,
      y = "Trait Value"
    ) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8)
    )
  if (.perc) {
    return (p + scale_y_continuous(
      labels = scales::percent_format(),
      expand = expansion(mult = c(0.01, 0.08))
    ))
  } else {
    return (p + scale_y_continuous(
      expand = expansion(mult = c(0.01, 0.08))
    ))
  }
}

plots <- plot_data %>%
  mutate(plot = list(plot_boxplot(data, trait, description, percent)))
p <- reduce(plots$plot, `+`) + plot_layout(nrow = 1, axes = "collect_y")
# tgutil::ggpreview(p, width = 12, height = 2.5)
ggsave("results/figures/diff_analysis/trait_boxplots.pdf", p, width = 12, height = 2.5)

write_csv(plot_data, "results/source_data/Figure_3k.csv")
