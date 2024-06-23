library(tidyverse)
library(ggrepel)
library(patchwork)

ancova_result <- read_csv(snakemake@input[[1]])

plot_dot <- function(data, .title) {
  data %>%
    mutate(label = if_else(p.adj < 0.05, glycan, "")) %>%
    ggplot(aes(glycan, -log10(p.adj))) +
    geom_point(aes(alpha = p.adj < 0.05), size = 2, color = "steelblue") +
    geom_text_repel(aes(label = label)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(x = "Derived trait", y = expression(paste("-log"[10], " p-value")), title = .title) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) +
    scale_y_continuous(
      limits = c(0, 10),
      expand = expansion(mult = c(0, 0.05))
    )
}

age_p <- ancova_result %>%
  filter(Effect == "age") %>%
  plot_dot("Effect: Age")

sex_p <- ancova_result %>%
  filter(Effect == "sex") %>%
  plot_dot("Effect: Sex")

p <- age_p | sex_p
# tgutil::ggpreview(width = 10, height = 5)
ggsave(snakemake@output[[1]], plot = p, width = 10, height = 5)