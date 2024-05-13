source("renv/activate.R")

library(tidyverse)
library(ggbeeswarm)
library(ggsignif)
library(rstatix)

# predictions <- read_csv("results/data/ml/predictions.csv")
# groups <- read_csv("results/data/prepared/groups.csv") |> 
#   mutate(group = factor(group, levels = c("H", "M", "Y", "C")))

predictions <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]]) |> 
  mutate(group = factor(group, levels = c("H", "M", "Y", "C")))

predictions <- predictions |> 
  inner_join(groups, by = "sample")

ggplot(predictions, aes(group, probability)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_quasirandom(aes(color = group)) +
  geom_signif(
    comparisons = list(c("H", "C"), c("M", "C"), c("Y", "C")),
    step_increase = 0.14,
    textsize = 3.5
  ) +
  annotate(
    "text", x = 0.5, y = 0.5, label = "Classification\nThreshold", 
    hjust = 0, size = 3, color = "grey50"
  ) +
  labs(x = "Group", y = "Predicted Probability") +
  theme_classic() +
  theme(legend.position = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_color_manual(values = c("H" = "#7A848D", "M" = "#A2AFA6", "Y" = "#FEC37D", "C" = "#CC5F5A"))
# tgutil::ggpreview(width = 3, height = 3)
ggsave(snakemake@output[[1]], width = 3, height = 3)