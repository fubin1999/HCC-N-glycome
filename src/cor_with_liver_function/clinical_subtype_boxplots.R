library(tidyverse)
library(rstatix)
library(ggsignif)
library(patchwork)

# traits <- read_csv("results/data/prepared/filtered_derived_traits.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")

traits <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clinical <- read_csv(snakemake@input[[3]])

traits <- traits %>%
  mutate(`CA3+CA4` = CA3 + CA4)

selected_traits <- c("TF", "CG", "CGS", "TB", "CA3+CA4")

trait_definitions <- tribble(
  ~trait, ~trait_definition, ~trait_field,
  "TF", "Prop. of Fucosylated Glycans", "Fucosylation",
  "CG", "Average No. of\nGalactoses per Antenna", "Galactosylation",
  "CGS", "Average No. of\nSia. per Gal.", "Sialylation",
  "TB", "Prop. of Bisecting Glycans", "Bisecting",
  "CA3+CA4", "Prop. of\nHigh-branching Glycans", "Branching"
)


subtypes <- clinical %>%
  mutate(
    sample = sample,
    child_pugh_subtype = if_else(child_pugh == "A", "CP A", "CP B/C"),
    ALBI_stage_subtype = if_else(ALBI_stage == "I", "ALBI I", "ALBI II/III"),
    .keep = "none"
  ) %>%
  pivot_longer(-sample, names_to = "subtype_by", values_to = "subtype", names_pattern = "(.*)_subtype") %>%
  filter(!is.na(subtype))

plot_data <- groups %>%
  filter(group != "QC") %>%
  left_join(subtypes, by = "sample") %>%
  left_join(traits %>% select(sample, all_of(selected_traits)), by = "sample") %>%
  pivot_longer(all_of(selected_traits), names_to = "trait", values_to = "trait_value")

draw_boxplot <- function (data, .ylab, .title) {
  ggplot(data, aes(subtype, trait_value)) +
    geom_boxplot(aes(color = subtype), outlier.alpha = 0) +
    geom_jitter(aes(color = subtype), width = 0.25, size = 1, shape = 16, alpha = 0.5) +
    geom_signif(
      comparisons = list(unique(data$subtype)),
      map_signif_level = function(p) sprintf("p = %.2g", p),
      size = 0.3, tip_length = 0, vjust = -0.2,
      test = "t.test"
    ) +
    scale_color_manual(values = c("#275D87", "#D26F32")) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.15)),
      labels = scales::percent_format()
    ) +
    labs(title = .title, y = .ylab) +
    guides(color = "none") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12)
    )
}

plot_df <- plot_data %>%
  nest_by(trait, subtype_by) %>%
  left_join(trait_definitions, by = "trait") %>%
  mutate(plot = list(draw_boxplot(data, trait_definition, trait_field))) %>%
  select(trait, subtype_by, plot) %>%
  arrange(subtype_by)

n_traits <- n_distinct(plot_df$trait)
n_subtypes <- n_distinct(plot_df$subtype_by)
plot_width <- n_traits * 2
plot_height <- n_subtypes * 2.5
p <- reduce(plot_df$plot, `+`) + plot_layout(nrow = n_subtypes, ncol = n_traits)

# tgutil::ggpreview(plot = p, width = plot_width, height = plot_height)
ggsave(snakemake@output[[1]], plot = p, width = plot_width, height = plot_height)