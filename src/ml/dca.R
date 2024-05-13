source("renv/activate.R")

library(tidyverse)

predictions <- read_csv("results/data/ml/predictions.csv")

thresholds <- seq(0, 0.95, by = 0.01)
prevalence <- mean(predictions$target)
net_benefit <- map_df(thresholds, function(threshold) {
  tp = sum(predictions$target == TRUE & predictions$probability >= threshold)
  fp = sum(predictions$target == FALSE & predictions$probability >= threshold)
  n = nrow(predictions)
  nb = (tp/n) - (fp/n) * (threshold / (1 - threshold))
  
  all_treat_nb = prevalence - (1 - prevalence) * (threshold / (1 - threshold))
  
  tibble(threshold = threshold, net_benefit = nb, all_treat = all_treat_nb, no_treat = 0)
})


ggplot(net_benefit, aes(x = threshold)) +
  geom_line(aes(y = net_benefit), color = "#CC5F5A") +
  geom_line(aes(y = all_treat), linetype = "dashed", color = "#7A848D") +
  geom_line(aes(y = no_treat), linetype = "dashed", color = "#A2AFA6") +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  labs(x = "Probability Threshold", y = "Net Benefit", title = "Decision Curve Analysis") +
  ylim(-0.1, 0.5) +
  theme(panel.grid = element_blank())

