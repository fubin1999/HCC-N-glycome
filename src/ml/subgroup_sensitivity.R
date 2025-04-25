library(tidyverse)


preds <- read_csv("results/data/ml/preds.csv")
clinical <- read_csv("results/data/prepared/clinical.csv")


AFP_neg_prepared <- preds %>% 
  mutate(model = case_match(
    model,
    "global" ~ "Model 4",
    "HCC_HC" ~ "Model 1",
    "HCC_CHB" ~ "Model 2",
    "HCC_LC" ~ "Model 3",
  )) |> 
  filter(dataset == "test", true) %>% 
  select(sample, pred, model) %>% 
  left_join(clinical %>% select(sample, AFP), by = "sample") %>% 
  filter(!is.na(AFP)) %>%
  mutate(
    AFP_neg_10 = AFP < 10,
    AFP_neg_20 = AFP < 20,
    AFP_neg_100 = AFP < 100,
    AFP_neg_200 = AFP < 200,
    AFP_neg_400 = AFP < 400,
    .keep = "unused"
  ) %>% 
  pivot_longer(
    cols = starts_with("AFP_neg"),
    names_to = "cutoff",
    values_to = "AFP_neg",
    names_prefix = "AFP_neg_"
  ) %>% 
  summarise(sensitivity = sum(pred & AFP_neg) / sum(AFP_neg), .by = c(model, cutoff)) %>% 
  mutate(
    cutoff = factor(as.integer(cutoff), levels = c(10, 20, 100, 200, 400)),
    model = factor(model, levels = c("Model 1", "Model 2", "Model 3", "Model 4"))
    )

ggplot(AFP_neg_prepared, aes(model, sensitivity, fill = cutoff)) +
  geom_col(position = "dodge") +
  labs(x = "Model", y = "Sensitivity", fill = "AFP Cutoff") +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
tgutil::ggpreview(width = 4, height = 3)
ggsave("results/figures/ml/AFP_neg_sensitivity.pdf", width = 4, height = 3)


TNM_prepared <- preds %>% 
  filter(dataset == "test", true) %>% 
  select(sample, pred, model) %>% 
  left_join(clinical %>% select(sample, TNM_stage), by = "sample") %>% 
  filter(!is.na(TNM_stage)) %>% 
  summarise(sensitivity = mean(pred), .by = c(model, TNM_stage)) %>% 
  mutate(
    model = factor(model, levels = c("HCC_HC", "HCC_CHB", "HCC_LC", "global")),
    TNM_stage = factor(TNM_stage, levels = c("I", "II", "III", "IV"))
  )

ggplot(TNM_prepared, aes(model, sensitivity, fill = TNM_stage)) +
  geom_col(position = "dodge") +
  labs(x = "Model", y = "Sensitivity", fill = "TNM Stage") +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
tgutil::ggpreview(width = 4, height = 3)
ggsave("results/figures/ml/TNM_sensitivity.pdf", width = 4, height = 3)
