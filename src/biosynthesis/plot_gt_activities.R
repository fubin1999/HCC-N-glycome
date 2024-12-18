library(tidyverse)


# activities <- read_csv("results/data/biosynthesis/GT_activities.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

activities <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])


prepared <- activities %>%
  pivot_longer(-glycotransferase, names_to = "sample", values_to = "activity") %>%
  left_join(groups, by = "sample") %>%
  filter(!group == "QC") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

p <- prepared %>%
  mutate(activity = scale(activity), .by = glycotransferase) %>%
  summarise(mean_activity = mean(activity), .by = c(group, glycotransferase)) %>%
  mutate(text_color = if_else(abs(mean_activity) > 0.6, "white", "black")) %>%
  ggplot(aes(group, glycotransferase, fill = mean_activity)) +
  geom_point(size = 5, color = "grey30", shape = 21) +
  labs(x = "Glycotransferase", y = "Group", fill = "Activity") +
  scale_fill_gradient2(high = "#D26F32", low = "#275D87", mid = "white") +
  theme_void() +
  theme(
    axis.text.y = element_text(),
    legend.position = "bottom",
  )

# tgutil::ggpreview(width = 1.5, height = 3)
ggsave(snakemake@output[[1]], p, width = 1.5, height = 3)