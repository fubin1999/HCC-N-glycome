source("renv/activate.R")

library(tidyverse)
library(ggbeeswarm)
library(ggsignif)
library(rstatix)

predictions <- read_csv("results/data/ml/predictions.csv")
groups <- read_csv("results/data/prepared/groups.csv") |>
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

predictions <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]]) |> 
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

predictions <- predictions |> 
  inner_join(groups, by = "sample")

ggplot(predictions, aes(group, probability)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_quasirandom(aes(color = group)) +
  geom_signif(
    comparisons = list(c("HC", "HCC"), c("CHB", "HCC"), c("LC", "HCC")),
    step_increase = 0.14,
    textsize = 3.5
  ) +
  annotate(
    "text", x = 0.5, y = 0.5, label = "Classification\nThreshold", 
    hjust = 0, size = 3, color = "grey50"
  ) +
  labs(x = "Group", y = "Predicted Probability") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_color_manual(values = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A"))
# tgutil::ggpreview(width = 3, height = 3)
ggsave(snakemake@output[[1]], width = 3, height = 3)