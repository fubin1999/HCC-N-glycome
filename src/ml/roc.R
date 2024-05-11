source("renv/activate.R")

library(tidyverse)
library(pROC)
library(patchwork)

# Read data-----
# predictions <- read_csv("results/data/ml/predictions.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# AFP_data <- read_csv("results/data/prepared/clinical.csv") |> 
#   select(sample, AFP)
# other_glycan_markers <- read_csv("results/data/prepared/other_glycan_markers.csv")

predictions <- read_csv(snakemake@input[["predictions"]])
groups <- read_csv(snakemake@input[["groups"]])
AFP_data <- read_csv(snakemake@input[["clinical"]]) |> 
  select(sample, AFP)
other_glycan_markers <- read_csv(snakemake@input[["other_glycan_markers"]])

predictions <- predictions |> 
  left_join(groups, by = "sample") |> 
  left_join(AFP_data, by = "sample") |> 
  left_join(other_glycan_markers, by = "sample")

# Draw ROC curve-----
draw_roc_curve <- function(glycan_roc, AFP_roc, Gtest_roc) {
  ggroc(list(
    `HCC Fusion` = glycan_roc, 
    AFP = AFP_roc,
    `G-test` = Gtest_roc
  )) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey") +
    labs(color = "") +
    theme_bw() +
    theme(
      legend.position = c(0.65, 0.2),
      panel.grid = element_blank(),
      legend.title = element_blank()
    ) +
    scale_color_manual(values = c("#CC5F5A", "#7A848D", "#A2AFA6"))
}

global_glycan_roc <- roc(predictions, response = "target", predictor = "probability", ci = TRUE)
global_AFP_roc <- roc(predictions, response = "target", predictor = "AFP", ci = TRUE)
global_Gtest_roc <- roc(predictions, response = "target", predictor = "Gtest", ci = TRUE)
global_p <- draw_roc_curve(global_glycan_roc, global_AFP_roc, global_Gtest_roc)

HC_glycan_roc <- roc(predictions |> filter(group %in% c("H", "C")), response = "target", predictor = "probability", ci = TRUE)
HC_AFP_roc <- roc(predictions |> filter(group %in% c("H", "C")), response = "target", predictor = "AFP", ci = TRUE)
HC_Gtest_roc <- roc(predictions |> filter(group %in% c("H", "C")), response = "target", predictor = "Gtest", ci = TRUE)
HC_p <- draw_roc_curve(HC_glycan_roc, HC_AFP_roc, HC_Gtest_roc)

MC_glycan_roc <- roc(predictions |> filter(group %in% c("M", "C")), response = "target", predictor = "probability", ci = TRUE)
MC_AFP_roc <- roc(predictions |> filter(group %in% c("M", "C")), response = "target", predictor = "AFP", ci = TRUE)
MC_Gtest_roc <- roc(predictions |> filter(group %in% c("M", "C")), response = "target", predictor = "Gtest", ci = TRUE)
MC_p <- draw_roc_curve(MC_glycan_roc, MC_AFP_roc, MC_Gtest_roc)

YC_glycan_roc <- roc(predictions |> filter(group %in% c("Y", "C")), response = "target", predictor = "probability", ci = TRUE)
YC_AFP_roc <- roc(predictions |> filter(group %in% c("Y", "C")), response = "target", predictor = "AFP", ci = TRUE)
YC_Gtest_roc <- roc(predictions |> filter(group %in% c("Y", "C")), response = "target", predictor = "Gtest", ci = TRUE)
YC_p <- draw_roc_curve(YC_glycan_roc, YC_AFP_roc, YC_Gtest_roc)

final_p <- global_p | HC_p | MC_p | YC_p
ggsave(snakemake@output[[1]], final_p, width = 16, height = 4)

# Save results-----
auc_df <- tribble(
  ~between, ~predictor, ~AUC, ~ci_lower, ~ci_upper,
  "Global", "HCC Fusion", global_glycan_roc$auc[1], global_glycan_roc$ci[1], global_glycan_roc$ci[3],
  "Global", "AFP", global_AFP_roc$auc[1], global_AFP_roc$ci[1], global_AFP_roc$ci[3],
  "Global", "G-test", global_Gtest_roc$auc[1], global_Gtest_roc$ci[1], global_Gtest_roc$ci[3],
  "HvsC", "HCC Fusion", HC_glycan_roc$auc[1], HC_glycan_roc$ci[1], HC_glycan_roc$ci[3],
  "HvsC", "AFP", HC_AFP_roc$auc[1], HC_AFP_roc$ci[1], HC_AFP_roc$ci[3],
  "HvsC", "G-test", HC_Gtest_roc$auc[1], HC_Gtest_roc$ci[1], HC_Gtest_roc$ci[3],
  "MvsC", "HCC Fusion", MC_glycan_roc$auc[1], MC_glycan_roc$ci[1], MC_glycan_roc$ci[3],
  "MvsC", "AFP", MC_AFP_roc$auc[1], MC_AFP_roc$ci[1], MC_AFP_roc$ci[3],
  "MvsC", "G-test", MC_Gtest_roc$auc[1], MC_Gtest_roc$ci[1], MC_Gtest_roc$ci[3],
  "YvsC", "HCC Fusion", YC_glycan_roc$auc[1], YC_glycan_roc$ci[1], YC_glycan_roc$ci[3],
  "YvsC", "AFP", YC_AFP_roc$auc[1], YC_AFP_roc$ci[1], YC_AFP_roc$ci[3],
  "YvsC", "G-test", YC_Gtest_roc$auc[1], YC_Gtest_roc$ci[1], YC_Gtest_roc$ci[3]
)
write_csv(auc_df, snakemake@output[[2]])
