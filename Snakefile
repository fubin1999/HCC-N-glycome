PREPARED_DIR = "results/data/prepared/"
RAW_ABUNDANCE = PREPARED_DIR + "raw_abundance.csv"
PROCESSED_ABUNDANCE = PREPARED_DIR + "processed_abundance.csv"
ALL_DERIVED_TRAITS = PREPARED_DIR + "derived_traits.csv"
FILTERED_DERIVED_TRAITS = PREPARED_DIR + "filtered_derived_traits.csv"
META_PROPERTIES = PREPARED_DIR + "meta_properties.csv"
GROUPS = PREPARED_DIR + "groups.csv"
CLINICAL = PREPARED_DIR + "clinical.csv"

rule all:
    input:
        # ===== Others =====
        "results/data/prepared/raw_abundance_full.csv",
        "results/data/clinical/AFP_cutoff.csv",
        "results/figures/clinical/AFP_cutoff.pdf",
        "results/figures/clinical/clinical_heatmap.pdf",
        "results/figures/SNFG/complete_svg/",
        "results/figures/SNFG/compact_svg/",
        "results/figures/SNFG/complete_pdf/",
        "results/figures/SNFG/compact_pdf/",

        # ===== Data Quality Figures =====
        "results/figures/data_quality/batch_effect_pca.pdf",
        "results/figures/data_quality/glycan_count_venn_all.pdf",
        "results/figures/data_quality/glycan_count_venn_with_structures.pdf",
        "results/figures/data_quality/glycan_count_venn_confident.pdf",
        "results/figures/data_quality/glycan_property_heatmap.pdf",
        "results/figures/data_quality/glycan_property_barplots.pdf",
        "results/figures/data_quality/glycan_property_boxplots.pdf",
        "results/figures/data_quality/QC_corrplot.pdf",

        # ===== Differential Analysis Data =====
        "results/data/diff_analysis/detect_rate_diff.csv",
        "results/data/diff_analysis/glycan_ancova.csv",
        "results/data/diff_analysis/glycan_post_hoc.csv",
        "results/data/diff_analysis/glycan_fold_change.csv",
        "results/data/diff_analysis/trait_ancova.csv",
        "results/data/diff_analysis/trait_post_hoc.csv",
        "results/data/diff_analysis/trait_fold_change.csv",

        # ===== Differential Analysis Figures =====
        "results/figures/diff_analysis/detect_rate_diff.pdf",
        "results/figures/diff_analysis/glycan_diff_rose_plot.pdf",
        "results/figures/diff_analysis/glycan_violin_plots.pdf",
        "results/figures/diff_analysis/glycan_diff_bubble.pdf",
        "results/figures/diff_analysis/glycan_heatmap.pdf",
        "results/figures/diff_analysis/glycan_volcanos.pdf",
        "results/figures/diff_analysis/glycan_compare_FC.pdf",
        "results/figures/diff_analysis/glycan_compare_p.pdf",
        "results/figures/diff_analysis/glycan_confounders.pdf",
        "results/figures/diff_analysis/glycan_pca.pdf",
        "results/figures/diff_analysis/trait_confounders.pdf",
        "results/figures/diff_analysis/trait_boxplots/",
        "results/figures/diff_analysis/trait_diff_bubble.pdf",
        "results/figures/diff_analysis/trait_heatmap.pdf",

        # ===== Glycan Coexpression Module Data =====
        "results/data/glycan_coexpr/glycan_clusters.csv",
        "results/data/glycan_coexpr/eigen_glycans.csv",
        "results/data/glycan_coexpr/cluster_ancova.csv",
        "results/data/glycan_coexpr/cluster_post_hoc.csv",
        "results/data/glycan_coexpr/cluster_cor_with_clinical.csv",

        # ===== Glycan Coexpression Module Figures =====
        "results/figures/glycan_coexpr/cluster_glycan_heatmap.pdf",
        "results/figures/glycan_coexpr/glycan_cluster_trends.pdf",
        "results/figures/glycan_coexpr/glycan_property_heatmap.pdf",
        "results/figures/glycan_coexpr/cluster_corrplot.pdf",
        "results/figures/glycan_coexpr/cluster_cor_with_clinical.pdf",
        "results/figures/glycan_coexpr/cor_inter_intra_GCM.pdf",

        # ===== Correlation with Clinical Data =====
        "results/data/cor_with_clinical/glycan_cor_with_AFP.csv",
        "results/data/cor_with_clinical/glycan_cor_with_clinical.csv",
        "results/data/cor_with_clinical/glycan_cor_with_clinical_per_group.csv",
        "results/data/cor_with_clinical/trait_cor_with_clinical.csv",
        "results/data/cor_with_clinical/trait_cor_with_clinical_per_group.csv",
        "results/data/cor_with_clinical/liver_function_model_r2.csv",
        "results/data/cor_with_clinical/liver_function_model_pred.csv",

        # ===== Correlation with Clinical Figures =====
        "results/figures/cor_with_clinical/glycan_cor_with_AFP_1.pdf",
        "results/figures/cor_with_clinical/glycan_cor_with_AFP_2.pdf",
        "results/figures/cor_with_clinical/trait_cor_with_clinical_all.pdf",
        "results/figures/cor_with_clinical/trait_cor_with_clinical_HC.pdf",
        "results/figures/cor_with_clinical/trait_cor_with_clinical_CHB.pdf",
        "results/figures/cor_with_clinical/trait_cor_with_clinical_LC.pdf",
        "results/figures/cor_with_clinical/trait_cor_with_clinical_HCC.pdf",
        "results/figures/cor_with_clinical/trait_clinical_subtype_boxplots.pdf",
        "results/figures/cor_with_clinical/glycan_cor_with_clinical.pdf",
        "results/figures/cor_with_clinical/liver_function_model_r2_venns.pdf",

        # ===== Motif Data =====
        # "results/data/GlyCompare_results/",
        # "results/data/motifs/motifs.csv",
        # "results/data/motifs/motif_structures.csv",
        # "results/data/motifs/ancova_result.csv",
        # "results/data/motifs/post_hoc_result.csv",

        # ===== Motif Figures =====
        # "results/figures/SNFG/glycomotifs/",

        # ===== TCGA Data =====
        "results/data/TCGA/dea_results.csv",
        "results/data/TCGA/consensus_cluster_result.csv",
        "results/data/TCGA/cluster_dea_results.csv",
        "results/data/TCGA/cox.csv",

        # ===== TCGA Figures =====
        "results/figures/TCGA/volcano.pdf",
        "results/figures/TCGA/MA_plot.pdf",
        "results/figures/TCGA/heatmap.pdf",
        "results/figures/TCGA/consensus_cluster/",
        "results/figures/TCGA/cluster_dea_heatmap.pdf",
        "results/figures/TCGA/clinical_heatmap.pdf",
        "results/figures/TCGA/cluster_KM.pdf",
        "results/figures/TCGA/single_gene_KM/",

        # ===== Machine Learning Data =====
        "results/data/ml/model_comparison.csv",
        "results/data/ml/mrmr_result.csv",
        "results/data/ml/predictions.csv",
        "results/data/ml/model_performance.csv",

        # ===== Machine Learning Figures =====
        "results/figures/ml/model_comparison_heatmap.pdf",
        "results/figures/ml/mrmr_auc.pdf",
        "results/figures/ml/roc_curves.pdf",
        "results/figures/ml/pr_curves.pdf",
        "results/figures/ml/calibration_curve.pdf",
        "results/figures/ml/confusion_matrix.pdf",
        "results/figures/ml/probability_boxplots.pdf",
        "results/figures/ml/complex_model_metrics_table.pdf",
        "results/figures/ml/simple_model_metrics_table.pdf",
        "results/figures/ml/forest_plot.pdf",
        "results/figures/ml/shap_summary.pdf",
        "results/figures/ml/shap_waterfall/",
        "results/figures/ml/decision_tree.svg"


# ==================== Prepare Data ====================
rule build_db:
    # Convert the serum glycan CSV file into a byonic database for GlyHunter.
    input:
        "data/human_serum_glycans.csv"
    output:
        "results/data/glycan_db.byonic"
    script:
        "src/prepare_data/make_db.py"

rule run_glyhunter:
    # Run GlyHunter on each mass list file.
    input:
        mass_list="data/mass_list/plate{no}.xlsx",
        db="results/data/glycan_db.byonic"
    output:
        directory("results/data/glyhunter_results/plate{no}/")
    shell:
        "glyhunter run {input.mass_list} -d {input.db} -o {output}"

rule run_glyhunter_default_db:
    # Run GlyHunter with default serum N-glycan database.
    input:
        "data/mass_list/plate{no}.xlsx"
    output:
        directory("results/data/glyhunter_results_full/plate{no}/")
    shell:
        "glyhunter run {input} -o {output}"

rule assign_maldi_pos:
    # Match the MALDI positions of GlyHunter results to samples.
    input:
        "data/MALDI_positions.csv",
        "results/data/glyhunter_results/plate{no}/"
    output:
        "results/data/data_per_plate/plate{no}.csv"
    script:
        "src/prepare_data/assign_maldi_pos.R"

rule assign_maldi_pos_full:
    # Match the MALDI positions of GlyHunter results to samples.
    input:
        "data/MALDI_positions.csv",
        "results/data/glyhunter_results_full/plate{no}/"
    output:
        "results/data/data_per_plate_full/plate{no}.csv"
    script:
        "src/prepare_data/assign_maldi_pos.R"

rule combine_data_per_plate:
    # Combine all GlyHunter results into a single file.
    input:
        expand("results/data/data_per_plate/plate{no}.csv", no=range(1, 9))
    output:
        RAW_ABUNDANCE
    script:
        "src/prepare_data/combine_data_per_plate.R"

rule combine_data_per_plate_full:
    # Combine all GlyHunter results into a single file.
    input:
        expand("results/data/data_per_plate_full/plate{no}.csv", no=range(1, 9))
    output:
        "results/data/prepared/raw_abundance_full.csv"
    script:
        "src/prepare_data/combine_data_per_plate.R"

rule prepare_groups:
    # Prepare the groups.
    input:
        plates="data/plates.csv"
    output:
        PREPARED_DIR + "unfiltered_groups.csv"
    script:
        "src/prepare_data/prepare_groups.R"

rule prepare_clinical:
    # Prepare the clinical information.
    input:
        clinical="data/clinical.csv",
        plates="data/plates.csv"
    output:
        PREPARED_DIR + "unfiltered_clinical.csv"
    script:
        "src/prepare_data/prepare_clinical.R"

rule preprocess:
    # Filter glycan, impute missing values, and normalize.
    input:
        RAW_ABUNDANCE,
        PREPARED_DIR + "unfiltered_groups.csv",
        PREPARED_DIR + "unfiltered_clinical.csv"
    output:
        PROCESSED_ABUNDANCE,
        GROUPS,
        CLINICAL
    script:
        "src/prepare_data/preprocess.R"

rule calculate_derived_traits:
    # Calculate derived traits using GlyTrait
    input:
        PROCESSED_ABUNDANCE,
        "data/human_serum_glycans.csv",
        "src/prepare_data/struc_builtin_formulas.txt"
    output:
        filtered_traits=FILTERED_DERIVED_TRAITS,
        all_traits=ALL_DERIVED_TRAITS,
        mp_table=META_PROPERTIES
    script:
        "src/prepare_data/derive_traits.py"


# ==================== Data Quality ====================
rule batch_effect_pca:
    # Draw PCA plot to check batch effect.
    input:
        PROCESSED_ABUNDANCE,
        "data/plates.csv"
    output:
        "results/figures/data_quality/batch_effect_pca.pdf"
    script:
        "src/data_quality/batch_effect.R"

rule glycan_count_venn:
    # Draw venn diagrams for glycan count per group.
    input:
        "results/data/prepared/raw_abundance_full.csv",
        "results/data/prepared/raw_abundance.csv",
        "results/data/prepared/processed_abundance.csv",
        "results/data/prepared/groups.csv"
    output:
        "results/figures/data_quality/glycan_count_venn_all.pdf",
        "results/figures/data_quality/glycan_count_venn_with_structures.pdf",
        "results/figures/data_quality/glycan_count_venn_confident.pdf"
    script:
        "src/data_quality/glycan_count_venn.R"

rule glycan_property_heatmap:
    # Draw heatmap to show all glycans' properties.
    input:
        PROCESSED_ABUNDANCE,
        RAW_ABUNDANCE,
        META_PROPERTIES
    output:
        "results/figures/data_quality/glycan_property_heatmap.pdf"
    script:
        "src/data_quality/glycan_property_heatmap.R"

rule glycan_property_detailed:
    # Draw barplots and boxplots for more detailed properties.
    # Barplots: for number of glycans in each category.
    # Boxplots: for relative abundance of glycans in each category.
    input:
        PROCESSED_ABUNDANCE,
        META_PROPERTIES
    output:
        "results/figures/data_quality/glycan_property_barplots.pdf",
        "results/figures/data_quality/glycan_property_boxplots.pdf"
    script:
        "src/data_quality/glycan_property_detailed.R"

rule QC_corrplot:
    # Draw correlation plot for QC samples.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS
    output:
        "results/figures/data_quality/QC_corrplot.pdf"
    script:
        "src/data_quality/QC_corrplot.R"


# ==================== Clinical Information ====================
rule AFP_cutoff:
    # Calculate classification metrics of AFP with different thresholds.
    input:
        CLINICAL,
        GROUPS
    output:
        "results/data/clinical/AFP_cutoff.csv",
        "results/figures/clinical/AFP_cutoff.pdf"
    script:
        "src/clinical/AFP_cutoff.R"

rule clinical_heatmap:
    # Draw heatmap to show clinical information.
    input:
        CLINICAL,
        GROUPS
    output:
        "results/figures/clinical/clinical_heatmap.pdf"
    script:
        "src/clinical/clinical_heatmap.R"


# ==================== Differential Analysis ====================
rule detection_rate_diff:
    # Perform Fisher's exact test to test detection rate difference of glycans.
    input:
        RAW_ABUNDANCE,
        GROUPS
    output:
        "results/data/diff_analysis/detect_rate_diff.csv",
        "results/figures/diff_analysis/detect_rate_diff.pdf"
    script:
        "src/diff_analysis/glycan_detect_rate.R"

rule glycan_ancova:
    # Perform ANCOVA for each glycan.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        clinical=CLINICAL
    output:
        ancova="results/data/diff_analysis/glycan_ancova.csv",
        posthoc="results/data/diff_analysis/glycan_post_hoc.csv"
    script:
        "src/diff_analysis/glycan_ancova.R"

rule glycan_fold_change:
    # Calculate fold changes for each glycan.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS
    output:
        "results/data/diff_analysis/glycan_fold_change.csv"
    script:
        "src/diff_analysis/glycan_fold_change.R"

rule glycan_diff_rose_plot:
    # Draw rose plot for differential glycans between each group pair.
    input:
        "results/data/diff_analysis/glycan_ancova.csv",
        "results/data/diff_analysis/glycan_post_hoc.csv"
    output:
        "results/figures/diff_analysis/glycan_diff_rose_plot.pdf"
    script:
        "src/diff_analysis/glycan_diff_rose_plot.R"

rule glycan_diff_bubble:
    # Draw bubble plot for p-values and fold changes of glycans.
    input:
        "results/data/diff_analysis/glycan_ancova.csv",
        "results/data/diff_analysis/glycan_post_hoc.csv",
        "results/data/diff_analysis/glycan_fold_change.csv",
    output:
        "results/figures/diff_analysis/glycan_diff_bubble.pdf"
    script:
        "src/diff_analysis/glycan_diff_bubble.R"

rule glycan_heatmap:
    # Draw mean heatmap for differential glycans.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS,
        "results/data/diff_analysis/glycan_ancova.csv"
    output:
        "results/figures/diff_analysis/glycan_heatmap.pdf"
    script:
        "src/diff_analysis/glycan_heatmap.R"

rule glycan_volcano:
    # Draw glycan volcano plots.
    input:
        "results/data/diff_analysis/glycan_ancova.csv",
        "results/data/diff_analysis/glycan_post_hoc.csv",
        "results/data/diff_analysis/glycan_fold_change.csv"
    output:
        "results/figures/diff_analysis/glycan_volcanos.pdf"
    script:
        "src/diff_analysis/glycan_volcano.R"

rule glycan_violin_plots:
    # Plot violin plots for all significant glycans.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        ancova_result="results/data/diff_analysis/glycan_ancova.csv"
    output:
        "results/figures/diff_analysis/glycan_violin_plots.pdf"
    script:
        "src/diff_analysis/glycan_violin_plots.R"

rule compare_glycan_FC:
    # Compare the fold changes of glycans between each group pair.
    input:
        "results/data/diff_analysis/glycan_fold_change.csv"
    output:
        "results/figures/diff_analysis/glycan_compare_FC.pdf"
    script:
        "src/diff_analysis/compare_FC.R"

rule compare_glycan_p:
    # Compare the p-values of glycans between each group pair.
    input:
        "results/data/diff_analysis/glycan_post_hoc.csv"
    output:
        "results/figures/diff_analysis/glycan_compare_p.pdf"
    script:
        "src/diff_analysis/compare_p.R"

rule glycan_confounders:
    # Plot dot plot for confounders' p-values.
    input:
        "results/data/diff_analysis/glycan_ancova.csv"
    output:
        "results/figures/diff_analysis/glycan_confounders.pdf"
    script:
        "src/diff_analysis/glycan_confounder_dot_plot.R"

rule glycan_pca:
    # Draw PCA plots for glycan abundance.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS
    output:
        "results/figures/diff_analysis/glycan_pca.pdf"
    script:
        "src/diff_analysis/glycan_pca.R"

rule trait_ancova:
    # Perform ANCOVA on derived traits.
    input:
        traits=FILTERED_DERIVED_TRAITS,
        groups=GROUPS,
        clinical=CLINICAL
    output:
        "results/data/diff_analysis/trait_ancova.csv",
        "results/data/diff_analysis/trait_post_hoc.csv"
    script:
        "src/diff_analysis/trait_ancova.R"

rule trait_fold_change:
    # Calculate fold changes for each derived trait.
    input:
        FILTERED_DERIVED_TRAITS,
        GROUPS
    output:
        "results/data/diff_analysis/trait_fold_change.csv"
    script:
        "src/diff_analysis/trait_fold_change.R"

rule trait_confounders:
    # Plot dot plot for confounders' p-values.
    input:
        "results/data/diff_analysis/trait_ancova.csv"
    output:
        "results/figures/diff_analysis/trait_confounders.pdf"
    script:
        "src/diff_analysis/trait_confounder_dot_plot.R"

rule trait_boxplots:
    # Draw boxplots for selected derived traits.
    input:
        ALL_DERIVED_TRAITS,
        GROUPS,
        "results/data/diff_analysis/trait_post_hoc.csv"
    output:
        directory("results/figures/diff_analysis/trait_boxplots/")
    script:
        "src/diff_analysis/trait_boxplots.R"

rule trait_diff_bubble:
    # Draw bubble plot for differential derived traits.
    input:
        "results/data/diff_analysis/trait_post_hoc.csv"
    output:
        "results/figures/diff_analysis/trait_diff_bubble.pdf"
    script:
        "src/diff_analysis/trait_diff_bubble.R"

rule trait_heatmap:
    # Draw heatmap for derived traits.
    input:
        FILTERED_DERIVED_TRAITS,
        GROUPS,
        "results/data/diff_analysis/trait_post_hoc.csv"
    output:
        "results/figures/diff_analysis/trait_heatmap.pdf"
    script:
        "src/diff_analysis/trait_heatmap.R"
    

# ==================== Glycan Coexpression Module ====================
rule glycan_cluster_heatmap:
    # Draw heatmap for differential glycans.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        ancova_result="results/data/diff_analysis/glycan_ancova.csv",
        mp_table=META_PROPERTIES
    output:
        "results/figures/glycan_coexpr/cluster_glycan_heatmap.pdf",
        "results/data/glycan_coexpr/glycan_clusters.csv"
    script:
        "src/glycan_coexpr/cluster_heatmap.R"

rule glycan_cor_within_cluster:
    # Correlation analysis for glycans within the same GCM or with different GCMs.
    input:
        PROCESSED_ABUNDANCE,
        "results/data/glycan_coexpr/glycan_clusters.csv"
    output:
        "results/figures/glycan_coexpr/cor_inter_intra_GCM.pdf"
    script:
        "src/glycan_coexpr/glycan_cor_within_cluster.R"

rule eigen_glycans:
    # Calculate the eigen value for each glycan cluster.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS,
        "results/data/glycan_coexpr/glycan_clusters.csv"
    output:
        "results/data/glycan_coexpr/eigen_glycans.csv"
    script:
        "src/glycan_coexpr/eigen_glycans.R"

rule cluster_ancova:
    # ANCOVA for glycan clusters.
    input:
        "results/data/glycan_coexpr/eigen_glycans.csv",
        GROUPS,
        CLINICAL
    output:
        "results/data/glycan_coexpr/cluster_ancova.csv",
        "results/data/glycan_coexpr/cluster_post_hoc.csv"
    script:
        "src/glycan_coexpr/cluster_ancova.R"

rule plot_cluster_trends:
    # Plot the alteration trends of glycan clusters from the heatmap about.
    input:
        "results/data/glycan_coexpr/eigen_glycans.csv",
        GROUPS,
        "results/data/glycan_coexpr/cluster_post_hoc.csv"
    output:
        "results/figures/glycan_coexpr/glycan_cluster_trends.pdf"
    script:
        "src/glycan_coexpr/plot_cluster_trends.R"

rule plot_cluster_properties:
    # Plot heatmap showing glycans' properties.
    input:
        "results/data/glycan_coexpr/glycan_clusters.csv",
        META_PROPERTIES
    output:
        "results/figures/glycan_coexpr/glycan_property_heatmap.pdf"
    script:
        "src/glycan_coexpr/cluster_property_heatmap.R"

rule cluster_corrplot:
    # Draw cluster corrplot.
    input:
        "results/data/glycan_coexpr/eigen_glycans.csv"
    output:
        "results/figures/glycan_coexpr/cluster_corrplot.pdf"
    script:
        "src/glycan_coexpr/cluster_cor.R"

rule cluster_cor_with_clinical:
    # Calculate and plot correlation of glycan clusters with clinical data.
    input:
        "results/data/glycan_coexpr/eigen_glycans.csv",
        GROUPS,
        CLINICAL
    output:
        "results/data/glycan_coexpr/cluster_cor_with_clinical.csv",
        "results/figures/glycan_coexpr/cluster_cor_with_clinical.pdf"
    script:
        "src/glycan_coexpr/cluster_cor_with_clinical.R"


# ==================== Correlation with Clinical Data ====================
rule glycan_cor_with_AFP:
    # Draw corrplot for glycans and AFP.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS,
        CLINICAL
    output:
        "results/data/cor_with_clinical/glycan_cor_with_AFP.csv",
        # There are too many glycans to display in a row,
        # and corrplot object could not be jointed by cowplot,
        # so the plot is splited into two.
        "results/figures/cor_with_clinical/glycan_cor_with_AFP_1.pdf",
        "results/figures/cor_with_clinical/glycan_cor_with_AFP_2.pdf"
    script:
        "src/cor_with_clinical/glycan_cor_with_AFP.R"

rule calcu_cor_with_clinical:
    # Correlation of glycan abundance and derived traits with clinical information.
    input:
        glycans=PROCESSED_ABUNDANCE,
        traits=FILTERED_DERIVED_TRAITS,
        groups=GROUPS,
        clinical=CLINICAL
    output:
        "results/data/cor_with_clinical/glycan_cor_with_clinical.csv",
        "results/data/cor_with_clinical/glycan_cor_with_clinical_per_group.csv",
        "results/data/cor_with_clinical/trait_cor_with_clinical.csv",
        "results/data/cor_with_clinical/trait_cor_with_clinical_per_group.csv"
    script:
        "src/cor_with_clinical/cor_with_clinical.R"

rule trait_corrplot_with_clinical:
    # Draw corrplot for derived traits with clinical information.
    input:
        "results/data/cor_with_clinical/trait_cor_with_clinical.csv",
        "results/data/cor_with_clinical/trait_cor_with_clinical_per_group.csv"
    output:
        all="results/figures/cor_with_clinical/trait_cor_with_clinical_all.pdf",
        HC="results/figures/cor_with_clinical/trait_cor_with_clinical_HC.pdf",
        CHB="results/figures/cor_with_clinical/trait_cor_with_clinical_CHB.pdf",
        LC="results/figures/cor_with_clinical/trait_cor_with_clinical_LC.pdf",
        HCC="results/figures/cor_with_clinical/trait_cor_with_clinical_HCC.pdf"
    script:
        "src/cor_with_clinical/trait_corrplot_with_clinical.R"

rule clinical_subtype_boxplots:
    # Draw boxplots for differential traits in clinical subtypes.
    input:
        FILTERED_DERIVED_TRAITS,
        GROUPS,
        CLINICAL
    output:
        "results/figures/cor_with_clinical/trait_clinical_subtype_boxplots.pdf"
    script:
        "src/cor_with_clinical/clinical_subtype_boxplots.R"

rule glycan_corrplot_with_clinical:
    # Draw corrplot for glycan abundance with clinical information.
    input:
        "results/data/cor_with_clinical/glycan_cor_with_clinical.csv",
    output:
        "results/figures/cor_with_clinical/glycan_cor_with_clinical.pdf"
    script:
        "src/cor_with_clinical/glycan_corrplot_with_clinical.R"

rule liver_function_model:
    # Build models to predict liver function.
    input:
        FILTERED_DERIVED_TRAITS,
        GROUPS,
        CLINICAL
    output:
        "results/data/cor_with_clinical/liver_function_model_r2.csv",
        "results/data/cor_with_clinical/liver_function_model_pred.csv"
    notebook:
        "src/cor_with_clinical/liver_function_model.ipynb"

rule liver_function_model_r2_venns:
    # Draw venn diagrams for R2 values of liver function models.
    input:
        "results/data/cor_with_clinical/liver_function_model_r2.csv"
    output:
        "results/figures/cor_with_clinical/liver_function_model_r2_venns.pdf"
    script:
        "src/cor_with_clinical/r2_venn.R"


# ==================== GlyCompare ====================
# GLYCOMPARE_CLI = "/Users/fubin/Python/glyCompareCT/glyCompareCT"
#
# rule prepare_for_glycompare:
#     # Prepare data format for GlyCompareCT.
#     input:
#         "data/glycan_structure_guess_linkage.csv"
#     output:
#         temp("results/data/glycompare_structures.csv")
#     run:
#         import pandas as pd
#         df = pd.read_csv(input[0])
#         df = df.rename(columns={"composition": "Name", "structure": "Glycan Structure"})
#         df.to_csv(output[0], index=False)
#
# rule run_glycompare:
#     # Run GlyCompareCT.
#     input:
#         PROCESSED_ABUNDANCE,
#         "results/data/glycompare_structures.csv"
#     output:
#         directory("results/data/GlyCompare_results/")
#     shell:
#         "{GLYCOMPARE_CLI} structure -a {input[0]} -v {input[1]} -o {output[0]} -p glycoCT -r N -c 8"
#
# rule plot_motif_SNFG:
#     # Draw SNFG cartoons for GlyCompare motifs.
#     input:
#         "results/data/GlyCompare/GlyCompare_output_data/GlyCompare_motif_annotation.csv"
#     output:
#         directory("results/figures/SNFG/glycomotifs/")
#     run:
#         import csv
#         from pathlib import Path
#         from glycowork.motif.draw import GlycoDraw
#
#         output_dir = output[0]
#         Path(output_dir).mkdir(exist_ok=True, parents=True)
#         with open(input[0], encoding='utf-8-sig') as fp:
#             reader = csv.reader(fp)
#             next(reader)
#             for motif, _, _, structure in reader:
#                 path = Path(output_dir) / f"{motif}.svg"
#                 GlycoDraw(structure, filepath=str(path), compact=True)
#
# rule tidy_glycompare_results:
#     # Convert the results of GlyCompare into tidy formats.
#     input:
#         "results/data/GlyCompare_results/GlyCompare_output_data/GlyCompare_motif_abd_table.csv",
#         "results/data/GlyCompare_results/GlyCompare_output_data/GlyCompare_motif_annotation.csv",
#     output:
#         "results/data/motifs/motifs.csv",
#         "results/data/motifs/motif_structures.csv"
#     run:
#         import pandas as pd
#
#         abund = pd.read_csv(input[0], index_col=0)
#         abund = abund.T
#         abund.index.name = "sample"
#         abund.to_csv(output[0])
#
#         struc = pd.read_csv(input[1], index_col=0)
#         struc = struc[["glycoCT"]]
#         struc = struc.rename(columns={"glycoCT": "structure"})
#         struc.index.name = "motif"
#         struc.to_csv(output[1])
#
# rule motif_ANCOVA:
#     # Perform ANCOVA on glycomotifs.
#     input:
#         "results/data/motifs/motifs.csv",
#         "results/data/prepared/groups.csv",
#         "results/data/prepared/clinical.csv",
#     output:
#         "results/data/motifs/ancova_result.csv",
#         "results/data/motifs/post_hoc_result.csv"
#     script:
#         "src/motif/ancova.R"


# ==================== TCGA Gene Expression ====================
rule download_TCGA:
    # Download TCGA gene expression data.
    output:
        "results/data/TCGA/query_summary.csv",
        directory("results/data/GDCdata"),
        "results/data/TCGA/prepared_data.rda"
    script:
        "src/TCGA/download.R"
        
rule DEA_TCGA:
    # Perform differential expression analysis on TCGA data.
    input:
        "results/data/TCGA/prepared_data.rda"
    output:
        "results/data/TCGA/dea_results.csv"
    script:
        "src/TCGA/DEA.R"

rule volcano_TCGA:
    # Draw volcano plot for TCGA gene expression.
    input:
        "results/data/TCGA/dea_results.csv",
        "data/glycogenes.csv"
    output:
        "results/figures/TCGA/volcano.pdf"
    script:
        "src/TCGA/volcano.R"

rule MA_TCGA:
    # Draw MA plot for TCGA gene expression.
    input:
        "results/data/TCGA/dea_results.csv",
        "data/glycogenes.csv"
    output:
        "results/figures/TCGA/MA_plot.pdf"
    script:
        "src/TCGA/MA_plot.R"

rule heatmap_TCGA:
    # Draw heatmap for differential glycogenes.
    input:
        "results/data/TCGA/prepared_data.rda",
        "results/data/TCGA/dea_results.csv",
        "data/glycogenes.csv"
    output:
        "results/figures/TCGA/heatmap.pdf"
    script:
        "src/TCGA/heatmap.R"

rule consensus_cluster:
    # Perform consensus clustering using glycogenes on tumor samples.
    input:
        "results/data/TCGA/prepared_data.rda",
        "data/glycogenes.csv"
    output:
        directory("results/figures/TCGA/consensus_cluster/"),
        "results/data/TCGA/consensus_cluster_result.csv"
    script:
        "src/TCGA/consensus_cluster.R"

rule DEA_cluster:
    # DEA on the two clusters identified by consensus clustering.
    input:
        "results/data/TCGA/prepared_data.rda",
        "results/data/TCGA/consensus_cluster_result.csv"
    output:
        "results/data/TCGA/cluster_dea_results.csv"
    script:
        "src/TCGA/DEA_cluster.R"

rule TCGA_cluster_heatmap:
    # Draw heatmap for different clusters.
    input:
        "results/data/TCGA/prepared_data.rda",
        "results/data/TCGA/consensus_cluster_result.csv",
        "data/glycogenes.csv"
    output:
        "results/figures/TCGA/cluster_dea_heatmap.pdf"
    script:
        "src/TCGA/heatmap_cluster.R"

rule TCGA_HCC_clinical_heatmap:
    # Draw clinical heatmap for TCGA HCC samples.
    input:
        "results/data/TCGA/prepared_data.rda",
        "results/data/TCGA/consensus_cluster_result.csv"
    output:
        "results/figures/TCGA/clinical_heatmap.pdf"
    script:
        "src/TCGA/clinical_heatmap_cluster.R"

rule TCGA_single_gene_survival:
    # Perform Cox PH on each glycogenes.
    input:
        "results/data/TCGA/prepared_data.rda",
        "data/glycogenes.csv"
    output:
        "results/data/TCGA/cox.csv",
        directory("results/figures/TCGA/single_gene_KM/")
    script:
        "src/TCGA/single_gene_survival.R"

rule TCGA_cluster_KM:
    # Draw KM Curve for two clusters.
    input:
        "results/data/TCGA/prepared_data.rda",
        "results/data/TCGA/consensus_cluster_result.csv"
    output:
        "results/figures/TCGA/cluster_KM.pdf"
    script:
        "src/TCGA/cluster_survival.R"


# ==================== Machine Learning ====================
rule prepare_data_for_ml:
    # Prepare data for ML and split the dataset.
    input:
        abundance=PROCESSED_ABUNDANCE,
        clinical=CLINICAL,
        groups=GROUPS
    output:
        train_data="results/data/ml/train_data.csv",
        test_data="results/data/ml/test_data.csv",
        feature_types="results/data/ml/feature_types.json"
    script:
        "src/ml/prepare_data.py"

rule compare_models:
    # Compare different machine learning models (including the HCC Fusion Classifier)
    # using 10-fold cross-validation.
    input:
        train_data="results/data/ml/train_data.csv",
        feature_types="results/data/ml/feature_types.json"
    output:
        "results/data/ml/model_comparison.csv"
    script:
        "src/ml/compare_models.py"

rule plot_compare_model_heatmap:
    # Draw the heatmap for model comparison
    input:
        "results/data/ml/model_comparison.csv"
    output:
        "results/figures/ml/model_comparison_heatmap.pdf"
    script:
        "src/ml/model_compare_heatmap.R"

rule mrmr_feature_selection:
    # Perform mRMR to select glycans.
    input:
        "results/data/ml/train_data.csv",
        "results/data/ml/feature_types.json"
    output:
        "results/data/ml/mrmr_result.csv"
    script:
        "src/ml/mrmr.py"

rule plot_mrmr:
    # Plot number of features vs. AUC for mRMR.
    input:
        "results/data/ml/mrmr_result.csv"
    output:
        "results/figures/ml/mrmr_auc.pdf"
    script:
        "src/ml/plot_mrmr.R"

rule make_predictions:
    # Predict on the test data and evalute the model.
    input:
        train_data="results/data/ml/train_data.csv",
        test_data="results/data/ml/test_data.csv",
        feature_types="results/data/ml/feature_types.json"
    output:
        "results/data/ml/predictions.csv"
    script:
        "src/ml/make_predictions.py"

rule model_roc_and_pr_curves:
    # Plot ROC curves and PR curves for two models and AFP.
    input:
        predictions="results/data/ml/predictions.csv",
        groups=GROUPS,
        clinical=CLINICAL
    output:
        "results/figures/ml/roc_curves.pdf",
        "results/figures/ml/pr_curves.pdf"
    script:
        "src/ml/pr_roc.R"

rule confusion_matrix:
    # Draw confusion matrix for the model.
    input:
        "results/data/ml/predictions.csv",
        "results/data/prepared/groups.csv"
    output:
        "results/figures/ml/confusion_matrix.pdf"
    script:
        "src/ml/confusion_matrix.R"

rule calibration_curve:
    # Draw the calibration curve of the model.
    input:
        "results/data/ml/predictions.csv"
    output:
        "results/figures/ml/calibration_curve.pdf"
    script:
        "src/ml/calibration_curve.R"

rule probability_boxplots:
    # Draw boxplots displaying probability distributions for each group.
    input:
        "results/data/ml/predictions.csv",
        "results/data/prepared/groups.csv"
    output:
        "results/figures/ml/probability_boxplots.pdf"
    script:
        "src/ml/score_boxplots.R"

rule evaluate_model:
    # Calculate the performance metrics of the model.
    input:
        "results/data/ml/predictions.csv",
        "results/data/prepared/groups.csv"
    output:
        "results/data/ml/model_performance.csv"
    script:
        "src/ml/metrics.R"

rule model_metrics_table:
    # Draw a table of model metrics.
    input:
        "results/data/ml/model_performance.csv"
    output:
        "results/figures/ml/complex_model_metrics_table.pdf",
        "results/figures/ml/simple_model_metrics_table.pdf"
    script:
        "src/ml/metrics_table.R"

rule forest_plot:
    # Draw a forest plot for the HCC Slim model.
    input:
        "results/data/ml/train_data.csv"
    output:
        "results/figures/ml/forest_plot.pdf"
    script:
        "src/ml/forest_plot.R"

rule shap:
    # SHAP on HCC Slim.
    input:
        "results/data/ml/train_data.csv",
        "results/data/ml/test_data.csv"
    output:
        "results/figures/ml/shap_summary.pdf",
        directory("results/figures/ml/shap_waterfall/")
    script:
        "src/ml/shap.py"

rule decision_tree:
    # Visualize a decision tree classifier.
    input:
        "results/data/ml/train_data.csv",
        "results/data/ml/test_data.csv"
    output:
        "results/figures/ml/decision_tree.svg"
    notebook:
        "src/ml/decision_tree.ipynb"


# ==================== Others ====================
rule draw_complete_glycans:
    # Draw glycans with glycowork with linkage information.
    input:
        "data/glycan_strucutre_best_resolution.csv"
    output:
        directory("results/figures/SNFG/complete_svg/")
    run:
        import csv
        from pathlib import Path
        from glycowork.motif.draw import GlycoDraw

        output_dir = output[0]
        Path(output_dir).mkdir(exist_ok=True, parents=True)
        with open(input[0], encoding='utf-8-sig') as fp:
            reader = csv.reader(fp)
            next(reader)
            for composition, structure in reader:
                path = Path(output_dir) / f"{composition}.svg"
                GlycoDraw(structure, filepath=str(path))

rule draw_compact_glycans:
    # Draw glycans with glycowork in the compact form.
    input:
        "data/glycan_structure_guess_linkage.csv"
    output:
        directory("results/figures/SNFG/compact_svg/")
    run:
        import csv
        from pathlib import Path
        from glycowork.motif.draw import GlycoDraw

        output_dir = output[0]
        Path(output_dir).mkdir(exist_ok=True, parents=True)
        with open(input[0], encoding='utf-8-sig') as fp:
            reader = csv.reader(fp)
            next(reader)
            for composition, structure in reader:
                path = Path(output_dir) / f"{composition}.svg"
                GlycoDraw(structure, filepath=str(path), compact=True)

rule convert_SNFG_svg_to_pdf:
    input:
        "results/figures/SNFG/{type}_svg/"
    output:
        directory("results/figures/SNFG/{type}_pdf/")
    run:
        from pathlib import Path
        from svglib.svglib import svg2rlg
        from reportlab.graphics import renderPDF
        output_path = Path(output[0])
        output_path.mkdir(parents=True)
        for svg_file in Path(input[0]).glob("*.svg"):
            name = svg_file.stem
            pdf_file = output_path / f"{name}.pdf"
            drawing = svg2rlg(str(svg_file))
            renderPDF.drawToFile(drawing, str(pdf_file))