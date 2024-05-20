import pandas as pd
from glytrait import Experiment, load_data


abundance_df = pd.read_csv(snakemake.input[0])
abundance_df = abundance_df.pivot(index="sample", columns="glycan", values="value")
abundance_df = abundance_df.reset_index()
abundance_df.columns.name = None
abundance_df = abundance_df.rename(columns={"sample": "Sample"})

structure_df = pd.read_csv(snakemake.input[1])
structure_df = pd.read_csv("data/human_serum_glycans.csv")
glycans = abundance_df.columns[1:].tolist()
structure_df = structure_df[structure_df["composition"].isin(glycans)]
structure_df = structure_df.reset_index(drop=True)
structure_df = structure_df.rename(columns={"composition": "GlycanID", "structure": "Structure"})

input_data = load_data(abundance_df, structure_df)
experiment = Experiment(input_data)
experiment.run_workflow(corr_threshold=0.9)

traits_df = experiment.filtered_derived_trait_table  # "Sample" as index, traits as columns
traits_df = traits_df.reset_index().melt(id_vars="Sample", var_name="trait", value_name="value")
traits_df = traits_df.rename(columns={"Sample": "sample"})
traits_df.to_csv(snakemake.output[0], index=False)

mp_df = experiment.meta_property_table
mp_df.index.name = "glycan"
mp_df.to_csv(snakemake.output[1], index=True)