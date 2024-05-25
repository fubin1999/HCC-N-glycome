PREPARED_DIR = "results/data/prepared/"
RAW_ABUNDANCE = PREPARED_DIR + "raw_abundance.csv"
PROCESSED_ABUNDANCE = PREPARED_DIR + "processed_abundance.csv"
GROUPS = PREPARED_DIR + "groups.csv"
CLINICAL = PREPARED_DIR + "clinical.csv"

rule all:
    input:
        # ===== Others =====
        "results/data/clinical/AFP_cutoff.csv",

        # ===== Data Quality Figures =====
        "results/figures/data_quality/batch_effect_pca.pdf",

        # ===== Differential Analysis Data =====
        "results/data/diff_analysis/ancova_for_glycans.csv",
        "results/data/diff_analysis/fold_change.csv",

        # ===== Differential Analysis Figures =====
        "results/figures/diff_analysis/diff_rose_plot.pdf",
        "results/figures/diff_analysis/diff_glycan_heatmap.pdf",
        "results/figures/diff_analysis/glycan_cluster_trends.pdf",
        "results/figures/diff_analysis/diff_bubble.pdf",
        "results/figures/diff_analysis/violin_plots.pdf",
        "results/figures/diff_analysis/diff_upset.pdf",

        # ===== Derived Traits Data =====
        "results/data/derived_traits/derived_traits.csv",
        "results/data/derived_traits/meta_properties.csv",

        # ===== Machine Learning Data =====
        "results/data/ml/model_comparison.csv",
        "results/data/ml/predictions.csv",
        "results/data/ml/roc_auc.csv",
        "results/data/ml/model_performance.csv",

        # ===== Machine Learning Figures =====
        "results/figures/ml/model_comparison_heatmap.pdf",
        "results/figures/ml/roc_curves.pdf",
        "results/figures/ml/calibration_curve.pdf",
        "results/figures/ml/confusion_matrix.pdf",
        "results/figures/ml/probability_boxplots.pdf",
        "results/figures/ml/model_metrics_table.pdf"


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

rule assign_maldi_pos:
    # Match the MALDI positions of GlyHunter results to samples.
    input:
        "data/MALDI_positions.csv",
        "results/data/glyhunter_results/plate{no}/"
    output:
        "results/data/data_per_plate/plate{no}.csv"
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

rule preprocess:
    # Filter glycan, impute missing values, and normalize.
    input:
        RAW_ABUNDANCE
    output:
        PROCESSED_ABUNDANCE
    script:
        "src/prepare_data/preprocess.R"

rule prepare_groups:
    # Prepare the groups.
    input:
        plates="data/plates.csv",
        abundance=PROCESSED_ABUNDANCE
    output:
        GROUPS
    script:
        "src/prepare_data/prepare_groups.R"

rule prepare_clinical:
    # Prepare the clinical information.
    input:
        clinical="data/clinical.csv",
        abundance=PROCESSED_ABUNDANCE,
        plates="data/plates.csv"
    output:
        CLINICAL
    script:
        "src/prepare_data/prepare_clinical.R"

rule derived_traits:
    # Calculate derived traits using GlyTrait
    input:
        PROCESSED_ABUNDANCE,
        "data/human_serum_glycans.csv"
    output:
        "results/data/derived_traits/derived_traits.csv",
        "results/data/derived_traits/meta_properties.csv"
    script:
        "src/derived_traits/derive_traits.py"

rule other_glycan_markers:
    # Calculate other HCC glycan markers
    input:
        PROCESSED_ABUNDANCE
    output:
        "results/data/prepared/other_glycan_markers.csv"
    script:
        "src/prepare_data/other_glycan_markers.R"


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


# ==================== Clinical Information ====================
rule AFP_cutoff:
    # Calculate classification metrics of AFP with different thresholds.
    input:
        CLINICAL,
        GROUPS
    output:
        "results/data/clinical/AFP_cutoff.csv"
    script:
        "src/clinical/AFP_cutoff.R"


# ==================== Differential Analysis ====================
rule ancova_for_glycans:
    # Perform ANCOVA for each glycan.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        clinical=CLINICAL
    output:
        ancova="results/data/diff_analysis/ancova_for_glycans.csv",
        posthoc="results/data/diff_analysis/posthoc_for_glycans.csv"
    script:
        "src/diff_analysis/ancova_for_glycans.R"

rule fold_change:
    # Calculate fold changes for each glycan.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS
    output:
        "results/data/diff_analysis/fold_change.csv"
    script:
        "src/diff_analysis/fold_change.R"

rule diff_rose_plot:
    # Draw rose plot for differential glycans between each group pair.
    input:
        "results/data/diff_analysis/posthoc_for_glycans.csv"
    output:
        "results/figures/diff_analysis/diff_rose_plot.pdf"
    script:
        "src/diff_analysis/rose_plot.R"

rule diff_glycan_heatmap:
    # Draw heatmap for differential glycans.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        ancova_result="results/data/diff_analysis/ancova_for_glycans.csv",
        mp_table="results/data/derived_traits/meta_properties.csv"
    output:
        "results/figures/diff_analysis/diff_glycan_heatmap.pdf",
        "results/data/diff_analysis/glycan_clusters.csv"
    script:
        "src/diff_analysis/heatmap.R"

rule diff_bubble:
    # Draw bubble plot for p-values and fold changes of glycans.
    input:
        post_hoc="results/data/diff_analysis/posthoc_for_glycans.csv",
        fold_change="results/data/diff_analysis/fold_change.csv",
        row_order="results/data/diff_analysis/glycan_clusters.csv"
    output:
        "results/figures/diff_analysis/diff_bubble.pdf"
    script:
        "src/diff_analysis/diff_bubble.R"

rule glycan_cluster_trends:
    # Plot the alteration trends of glycan clusters from the heatmap about.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        clusters="results/data/diff_analysis/glycan_clusters.csv"
    output:
        "results/figures/diff_analysis/glycan_cluster_trends.pdf"
    script:
        "src/diff_analysis/cluster_trends.R"

rule violin_plots:
    # Plot violin plots for all significant glycans.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        ancova_result="results/data/diff_analysis/ancova_for_glycans.csv"
    output:
        "results/figures/diff_analysis/violin_plots.pdf"
    script:
        "src/diff_analysis/violin_plots.R"

rule diff_upset:
    # Plot upset plot for the number of significant glycans
    # between each group paires.
    input:
        "results/data/diff_analysis/posthoc_for_glycans.csv"
    output:
        "results/figures/diff_analysis/diff_upset.pdf"
    script:
        "src/diff_analysis/diff_upset.R"


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

rule roc:
    # Calculate ROC AUCs for different groups and plot ROC curves.
    input:
        predictions="results/data/ml/predictions.csv",
        groups=GROUPS,
        clinical=CLINICAL,
        other_glycan_markers="results/data/prepared/other_glycan_markers.csv"
    output:
        "results/figures/ml/roc_curves.pdf",
        "results/data/ml/roc_auc.csv"
    script:
        "src/ml/roc.R"

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
        "results/figures/ml/model_metrics_table.pdf"
    script:
        "src/ml/metrics_table.R"