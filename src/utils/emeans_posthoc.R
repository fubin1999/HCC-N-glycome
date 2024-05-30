library(tidyverse)

post_hoc <- function (data, .formula, .target, .group_by) {
  n_samples <- data %>%
    distinct(sample, {{.target}}) %>%
    summarise(n = n(), .by = {{.target}})

  emmean_result <- data %>%
    group_by({{.group_by}}) %>%
    nest() %>%
    mutate(
      lm_model = map(data, ~ lm(.formula, data = .x)),
      emm = map(lm_model, ~ emmeans::emmeans(.x, ~ {{.target}})),
    ) %>%
    select({{.group_by}}, emm)

  post_hoc_result <- emmean_result %>%
    mutate(
      pairwise_comparison = map(emm, ~ emmeans::contrast(.x, "pairwise", adjust = "tukey")),
      result_df = map(pairwise_comparison, ~ as_tibble(.x))
    ) %>%
    select(-emm, -pairwise_comparison) %>%
    unnest(cols = result_df) %>%
    separate_wider_delim(contrast, delim = " - ", names = c("group1", "group2")) %>%
    rename(p.adj = p.value)

  emmean_values <- emmean_result %>%
    mutate(emm_summary = map(emm, ~ summary(.x))) %>%
    select({{.group_by}}, emm_summary) %>%
    unnest(cols = c(emm_summary)) %>%
    left_join(n_samples, by = join_by({{.target}})) %>%
    mutate(sd = SE * sqrt(n)) %>%
    select({{.group_by}}, {{.target}}, emmean, sd)

  effect_size <- post_hoc_result %>%
    select({{.group_by}}, group1, group2) %>%
    left_join(emmean_values, by = join_by({{.group_by}}, group1 == {{.target}})) %>%
    left_join(emmean_values, by = join_by({{.group_by}}, group2 == {{.target}}), suffix = c("1", "2")) %>%
    mutate(
      pooled_sd = sqrt((sd1^2 + sd2^2) / 2),
      cohens_d = (emmean2 - emmean1) / pooled_sd
    ) %>%
    select({{.group_by}}, group1, group2, cohens_d)

  post_hoc_result %>%
    left_join(effect_size, by = join_by({{.group_by}}, group1, group2))
}