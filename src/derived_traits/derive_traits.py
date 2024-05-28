import pandas as pd
from glytrait import Experiment, load_data
from glytrait.formula import load_formulas_from_file


abundance_df = pd.read_csv(snakemake.input[0])
abundance_df = abundance_df.rename(columns={"sample": "Sample"})

structure_df = pd.read_csv(snakemake.input[1])
# structure_df = pd.read_csv("data/human_serum_glycans.csv")
glycans = abundance_df.columns[1:].tolist()
structure_df = structure_df[structure_df["composition"].isin(glycans)]
structure_df = structure_df.reset_index(drop=True)
structure_df = structure_df.rename(columns={"composition": "GlycanID", "structure": "Structure"})

formula_file = snakemake.input[2]
formulas = load_formulas_from_file(formula_file)
input_data = load_data(abundance_df, structure_df)
experiment = Experiment(input_data=input_data)
experiment.run_workflow(formulas=formulas, corr_threshold=0.9)

filtered_traits_df = experiment.filtered_derived_trait_table  # "Sample" as index, traits as columns
filtered_traits_df = filtered_traits_df.reset_index()
filtered_traits_df = filtered_traits_df.rename(columns={"Sample": "sample"})
filtered_traits_df.to_csv(snakemake.output["filtered_traits"], index=False)

traits_df = experiment.derived_trait_table
traits_df = traits_df.reset_index()
traits_df = traits_df.rename(columns={"Sample": "sample"})
traits_df.to_csv(snakemake.output["all_traits"], index=False)

mp_df = experiment.meta_property_table
mp_df.index.name = "glycan"
mp_df.to_csv(snakemake.output["mp_table"], index=True)