library(tidyverse)
library(ggplate)

# Function-----
load_data = function(filename) {
  read_csv(filename) %>%
    extract(raw, into = "group", regex = "(.)", remove = F) %>%
    mutate(group = if_else(group == "Q", "QC", group))
}

draw_plate = function(data, version = "clr") {
  p = plate_plot(
    data = data,
    position = position,
    value = group,
    label = raw,
    label_size = 2.5,
    plate_size = 96,
    plate_type = "round"
  )
  if (version == "bw") {
    p +
      scale_fill_manual(values = rep("white", 5)) +
      theme(legend.position = "none", plot.title = element_blank())
  }
  else {
    p + 
      ggsci::scale_fill_npg() +
      theme(plot.title = element_blank())
  }
}


# Draw plates-----
library(patchwork)

data_5 = load_data("data/plates/plate7.csv")
data_6 = load_data("data/plates/plate8.csv")
p_5_bw = draw_plate(data_5, version = "bw")
p_6_bw = draw_plate(data_6, version = "bw")
p_5_bw + p_6_bw + plot_layout(ncol = 1)
tgutil::ggpreview(width = 140, height = 180, units = "mm")
ggsave(file="figures/plates/batch4.pdf", width = 140, height = 180, units = "mm")
