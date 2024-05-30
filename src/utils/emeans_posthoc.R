library(tidyverse)

post_hoc <- function (data, .formula, .target) {
  data %>%
    nest() %>%
    mutate(
      lm_model = map(data, ~ lm(.formula, data = .x)),
      emm = map(lm_model, ~ emmeans::emmeans(.x, ~ {{.target}})),
      pairwise_comparison = map(emm, ~ emmeans::contrast(.x, "pairwise", adjust = "tukey")),
      result_df = map(pairwise_comparison, ~ as_tibble(.x))
    ) %>%
    select(-data, -lm_model, -emm, -pairwise_comparison) %>%
    unnest(cols = result_df) %>%
    separate_wider_delim(contrast, delim = " - ", names = c("group1", "group2")) %>%
    rename(p.adj = p.value)
}