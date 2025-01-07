library(tidyverse)
library(easystats)
library(patchwork)


data <- read_csv("results/data/biosynthesis/GT_activities.csv")
groups <- read_csv("results/data/prepared/groups.csv")
clinical <- read_csv("results/data/prepared/clinical.csv")


prepared <- data %>%
  pivot_longer(-glycotransferase, names_to = "sample", values_to = "activity") %>%
  inner_join(groups, by = "sample") %>%
  inner_join(clinical %>% select(sample, sex, age, ALBI_score), by = "sample") %>%
  mutate(
    sex = factor(sex, levels = c("M", "F")),
    group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))
  )

models <- prepared %>%
  nest_by(glycotransferase) %>%
  mutate(model = list(lm(activity ~ sex + age + ALBI_score + group, data = data))) %>%
  select(-data)

model_performances <- models %>%
  mutate(performance = list(performance(model))) %>%
  select(-model) %>%
  unnest(performance)

model_coefficients <- models %>%
  mutate(coefficients = list(standardize_parameters(model, method = "refit", two_sd = TRUE))) %>%
  select(-model) %>%
  unnest(coefficients) %>%
  mutate(Parameter = case_match(
    Parameter,
    "sexF" ~ "Sex: F",
    "groupLC" ~ "Group: LC",
    "groupHCC" ~ "Group: HCC",
    "groupCHB" ~ "Group: CHB",
    "ALBI_score" ~ "ALBI Score",
    "age" ~ "Age",
    "(Intercept)" ~ "Intercept"
  )) %>% 
  mutate(Parameter = factor(
    Parameter, levels = rev(c("Group: HCC", "Group: LC", "Group: CHB", "ALBI Score",
                              "Sex: F", "Age", "Intercept"))
  )) %>% 
  ungroup()

plot_forest <- function(data) {
  ggplot(data, aes(Std_Coefficient, Parameter)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.1) +
    geom_point(aes(fill = Std_Coefficient), color = "black", size = 3, shape = 22) +
    scale_fill_gradient2(low = "#275D87", mid = "white", high = "#BE5800") +
    labs(x = "Std. Coefficient") +
    theme_classic() +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
}

plot_df <- model_coefficients %>%
  nest_by(glycotransferase) %>%
  mutate(plot = list(plot_forest(data) + ggtitle(glycotransferase)))
final_p <- reduce(plot_df$plot, `+`) + 
  plot_layout(nrow = 1, axes = "collect") &
  xlim(-1.1, 1.1)
tgutil::ggpreview(final_p, width = 15, height = 3)
ggsave("results/figures/biosynthesis/GT_activity_forest.pdf", final_p, width = 15, height = 3)
