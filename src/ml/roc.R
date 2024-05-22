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

# All data-----
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
      legend.position = c(0.65, 0.25),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.background = element_blank(),
      title = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#CC5F5A", "#7A848D", "#A2AFA6"))
}

global_glycan_roc <- roc(predictions, response = "target", predictor = "probability", ci = TRUE)
global_AFP_roc <- roc(predictions, response = "target", predictor = "AFP", ci = TRUE)
global_Gtest_roc <- roc(predictions, response = "target", predictor = "Gtest", ci = TRUE)
global_p <- draw_roc_curve(global_glycan_roc, global_AFP_roc, global_Gtest_roc) +
  ggtitle("Test Set: Control/HCC")

HC_glycan_roc <- roc(predictions |> filter(group %in% c("HC", "HCC")), response = "target", predictor = "probability", ci = TRUE)
HC_AFP_roc <- roc(predictions |> filter(group %in% c("HC", "HCC")), response = "target", predictor = "AFP", ci = TRUE)
HC_Gtest_roc <- roc(predictions |> filter(group %in% c("HC", "HCC")), response = "target", predictor = "Gtest", ci = TRUE)
HC_p <- draw_roc_curve(HC_glycan_roc, HC_AFP_roc, HC_Gtest_roc) +
  ggtitle("Test Set: HC/HCC")

MC_glycan_roc <- roc(predictions |> filter(group %in% c("CHB", "HCC")), response = "target", predictor = "probability", ci = TRUE)
MC_AFP_roc <- roc(predictions |> filter(group %in% c("CHB", "HCC")), response = "target", predictor = "AFP", ci = TRUE)
MC_Gtest_roc <- roc(predictions |> filter(group %in% c("CHB", "HCC")), response = "target", predictor = "Gtest", ci = TRUE)
MC_p <- draw_roc_curve(MC_glycan_roc, MC_AFP_roc, MC_Gtest_roc) +
  ggtitle("Test Set: CHB/HCC")

YC_glycan_roc <- roc(predictions |> filter(group %in% c("LC", "HCC")), response = "target", predictor = "probability", ci = TRUE)
YC_AFP_roc <- roc(predictions |> filter(group %in% c("LC", "HCC")), response = "target", predictor = "AFP", ci = TRUE)
YC_Gtest_roc <- roc(predictions |> filter(group %in% c("LC", "HCC")), response = "target", predictor = "Gtest", ci = TRUE)
YC_p <- draw_roc_curve(YC_glycan_roc, YC_AFP_roc, YC_Gtest_roc) +
  ggtitle("Test Set: LC/HCC")

# AFP negative samples-----
AFP_neg_samples <- predictions |> filter(AFP < 20)

draw_roc_curve <- function(glycan_roc, Gtest_roc) {
  ggroc(list(
    `HCC Fusion` = glycan_roc, 
    `G-test` = Gtest_roc
  )) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey") +
    labs(color = "") +
    theme_bw() +
    theme(
      legend.position = c(0.65, 0.2),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.background = element_blank(),
      title = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#CC5F5A", "#A2AFA6"))
}

AN_global_glycan_roc <- roc(AFP_neg_samples, response = "target", predictor = "probability", ci = TRUE)
AN_global_Gtest_roc <- roc(AFP_neg_samples, response = "target", predictor = "Gtest", ci = TRUE)
AN_global_p <- draw_roc_curve(AN_global_glycan_roc, AN_global_Gtest_roc) +
  ggtitle("Test Set: Control/HCC, AFP(-)")

AN_HC_glycan_roc <- roc(AFP_neg_samples |> filter(group %in% c("HC", "HCC")), response = "target", predictor = "probability", ci = TRUE)
AN_HC_Gtest_roc <- roc(AFP_neg_samples |> filter(group %in% c("HC", "HCC")), response = "target", predictor = "Gtest", ci = TRUE)
AN_HC_p <- draw_roc_curve(AN_HC_glycan_roc, AN_HC_Gtest_roc) +
  ggtitle("Test Set: HC/HCC, AFP(-)")

AN_MC_glycan_roc <- roc(AFP_neg_samples |> filter(group %in% c("CHB", "HCC")), response = "target", predictor = "probability", ci = TRUE)
AN_MC_Gtest_roc <- roc(AFP_neg_samples |> filter(group %in% c("CHB", "HCC")), response = "target", predictor = "Gtest", ci = TRUE)
AN_MC_p <- draw_roc_curve(AN_MC_glycan_roc, AN_MC_Gtest_roc) +
  ggtitle("Test Set: CHB/HCC, AFP(-)")

AN_YC_glycan_roc <- roc(AFP_neg_samples |> filter(group %in% c("LC", "HCC")), response = "target", predictor = "probability", ci = TRUE)
AN_YC_Gtest_roc <- roc(AFP_neg_samples |> filter(group %in% c("LC", "HCC")), response = "target", predictor = "Gtest", ci = TRUE)
AN_YC_p <- draw_roc_curve(AN_YC_glycan_roc, AN_YC_Gtest_roc) +
  ggtitle("Test Set: LC/HCC, AFP(-)")

final_p <- (global_p | HC_p | MC_p | YC_p) / (AN_global_p | AN_HC_p | AN_MC_p | AN_YC_p)
ggsave(snakemake@output[[1]], final_p, width = 12, height = 6)

# Save results-----
auc_df <- tribble(
  ~data, ~between, ~predictor, ~AUC, ~ci_lower, ~ci_upper,
  "all", "Control/HCC", "HCC Fusion", global_glycan_roc$auc[1], global_glycan_roc$ci[1], global_glycan_roc$ci[3],
  "all", "Control/HCC", "AFP", global_AFP_roc$auc[1], global_AFP_roc$ci[1], global_AFP_roc$ci[3],
  "all", "Control/HCC", "G-test", global_Gtest_roc$auc[1], global_Gtest_roc$ci[1], global_Gtest_roc$ci[3],
  "all", "HC/HCC", "HCC Fusion", HC_glycan_roc$auc[1], HC_glycan_roc$ci[1], HC_glycan_roc$ci[3],
  "all", "HC/HCC", "AFP", HC_AFP_roc$auc[1], HC_AFP_roc$ci[1], HC_AFP_roc$ci[3],
  "all", "HC/HCC", "G-test", HC_Gtest_roc$auc[1], HC_Gtest_roc$ci[1], HC_Gtest_roc$ci[3],
  "all", "CHB/HCC", "HCC Fusion", MC_glycan_roc$auc[1], MC_glycan_roc$ci[1], MC_glycan_roc$ci[3],
  "all", "CHB/HCC", "AFP", MC_AFP_roc$auc[1], MC_AFP_roc$ci[1], MC_AFP_roc$ci[3],
  "all", "CHB/HCC", "G-test", MC_Gtest_roc$auc[1], MC_Gtest_roc$ci[1], MC_Gtest_roc$ci[3],
  "all", "LC/HCC", "HCC Fusion", YC_glycan_roc$auc[1], YC_glycan_roc$ci[1], YC_glycan_roc$ci[3],
  "all", "LC/HCC", "AFP", YC_AFP_roc$auc[1], YC_AFP_roc$ci[1], YC_AFP_roc$ci[3],
  "all", "LC/HCC", "G-test", YC_Gtest_roc$auc[1], YC_Gtest_roc$ci[1], YC_Gtest_roc$ci[3],
  "AFP-", "Control/HCC", "HCC Fusion", AN_global_glycan_roc$auc[1], AN_global_glycan_roc$ci[1], AN_global_glycan_roc$ci[3],
  "AFP-", "Control/HCC", "G-test", AN_global_Gtest_roc$auc[1], AN_global_Gtest_roc$ci[1], AN_global_Gtest_roc$ci[3],
  "AFP-", "HC/HCC", "HCC Fusion", AN_HC_glycan_roc$auc[1], AN_HC_glycan_roc$ci[1], AN_HC_glycan_roc$ci[3],
  "AFP-", "HC/HCC", "G-test", AN_HC_Gtest_roc$auc[1], AN_HC_Gtest_roc$ci[1], AN_HC_Gtest_roc$ci[3],
  "AFP-", "CHB/HCC", "HCC Fusion", AN_MC_glycan_roc$auc[1], AN_MC_glycan_roc$ci[1], AN_MC_glycan_roc$ci[3],
  "AFP-", "CHB/HCC", "G-test", AN_MC_Gtest_roc$auc[1], AN_MC_Gtest_roc$ci[1], AN_MC_Gtest_roc$ci[3],
  "AFP-", "LC/HCC", "HCC Fusion", AN_YC_glycan_roc$auc[1], AN_YC_glycan_roc$ci[1], AN_YC_glycan_roc$ci[3],
  "AFP-", "LC/HCC", "G-test", AN_YC_Gtest_roc$auc[1], AN_YC_Gtest_roc$ci[1], AN_YC_Gtest_roc$ci[3]
)
write_csv(auc_df, snakemake@output[[2]])
