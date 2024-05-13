source("renv/activate.R")

library(tidyverse)
library(probably)

# predictions <- read_csv("results/data/ml/predictions.csv")
predictions <- read_csv(snakemake@input[[1]])

calibration_data <- predictions |> 
  mutate(prob_bin = cut_width(probability, width = 0.2, boundary = 0)) |>
  summarise(
    mean_prob = mean(probability), 
    event_rate = mean(target), 
    n = n(),
    lower_ci = binom.test(x = sum(target), n = n, conf.level = 0.9)$conf.int[1],
    upper_ci = binom.test(x = sum(target), n = n, conf.level = 0.9)$conf.int[2],
    .by = prob_bin
  )

ggplot(calibration_data, aes(mean_prob, event_rate)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.1, fill = "steelblue") +
  geom_point(color = "steelblue") +
  geom_line(color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggtitle("Calibration Curve") +
  labs(x = "Mean Predicted Probability", y = "Observed Event Rate") +
  lims(x = c(0, 1), y = c(0, 1)) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

# tgutil::ggpreview(width = 3, height = 3)
ggsave(snakemake@output[[1]], width = 3, height = 3)