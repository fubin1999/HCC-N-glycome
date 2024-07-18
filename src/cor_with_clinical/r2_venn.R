library(tidyverse)
library(ggVennDiagram)
library(patchwork)

# r2_result <- read_csv("results/data/cor_with_clinical/liver_function_model_r2.csv")
r2_result <- read_csv(snakemake@input[[1]])

extract_feature_table <- function (data, column, delim) {
  data %>%
    distinct({{column}}) %>%
    rename(features = {{column}}) %>%
    mutate(id = features) %>%
    tidyr::separate_longer_delim(features, delim = delim) %>%
    rename(feature = features, features = id) %>%
    mutate(used = TRUE) %>%
    complete(feature, features, fill = list(used = FALSE)) %>%
    pivot_wider(names_from = feature, values_from = used)
}

r2_feature_table <- extract_feature_table(r2_result, features, ",")

plot_r2_venn <- function (r2_data, .title) {
  # Generate an empty VennPlotData object
  data_list <- list()
  for (f in setdiff(colnames(r2_feature_table), "features")) {
    data_list[[f]] <- NA
  }
  venn_data <- process_data(Venn(data_list))

  # Use helper functions in ggVennDiagram to initialize the data.
  # This is necessary to get locations of the regions and labels (X, Y).
  set_edge_data <- venn_setedge(venn_data)
  set_label_data <- venn_setlabel(venn_data)
  region_label_data <- venn_regionlabel(venn_data)
  region_edge_data <- venn_regionedge(venn_data)

  # Update the data with the actual data
  venn_data_feature_table <- extract_feature_table(region_label_data, name, "/")
  set_label_data <- set_label_data %>%
    sjmisc::rec(name, rec = "CA3+CA4=Bra; CG=Gal; TB=Bis; TF=Fuc", suffix = "")
  region_label_data <- region_label_data %>%
    left_join(venn_data_feature_table, by = c("name" = "features")) %>%
    left_join(r2_data, by = setdiff(colnames(r2_feature_table), "features"))
  region_edge_data <- region_edge_data %>%
    left_join(venn_data_feature_table, by = c("name" = "features")) %>%
    left_join(r2_data, by = setdiff(colnames(r2_feature_table), "features"))

  # Plot the Venn diagram
  ggplot() +
    geom_polygon(aes(X, Y, fill = r2, group = id),
                 data = region_edge_data) +
    scale_fill_gradient(low = "white", high = "#D26F32", limits = c(0, 1)) +
    geom_path(aes(X, Y, group = id),
              color = "black",
              data = set_edge_data,
              show.legend = FALSE) +
    geom_text(aes(X, Y, label = name),
              data = set_label_data) +
    geom_text(aes(X, Y, label = scales::percent(r2, accuracy = 1)),
              data = region_label_data,
              size = 3
    ) +
    guides(fill = "none") +
    labs(caption = .title) +
    theme_void() +
    theme(
      plot.caption = element_text(size = 12, hjust = 0.5)
    )
}

plot_df <- r2_result %>%
  left_join(r2_feature_table, by = "features") %>%
  nest_by(target) %>%
  mutate(
    caption = str_glue("R² for predicting {target}"),
    plot = list(plot_r2_venn(data, caption))
  ) %>%
  select(target, plot)
p <- reduce(plot_df$plot, `+`) +
  plot_layout(nrow = 1)

# tgutil::ggpreview(width = 15, height = 3.5)
ggsave(snakemake@output[[1]], p, width = 15, height = 3)