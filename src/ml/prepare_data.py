import json

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split


# abundance = pd.read_csv("results/data/prepared/processed_abundance.csv")
# clinical = pd.read_csv("results/data/prepared/clinical.csv")
# groups = pd.read_csv("results/data/prepared/groups.csv")
abundance = pd.read_csv(snakemake.input["abundance"], index_col=0)
clinical = pd.read_csv(snakemake.input["clinical"])
groups = pd.read_csv(snakemake.input["groups"])

abundance = pd.DataFrame(np.log2(abundance.values), columns=abundance.columns, index=abundance.index)
clinical = clinical.drop(["sex", "age"], axis=1)
clinical = clinical.set_index("sample")
groups = groups.set_index("sample")
data = pd.merge(abundance, clinical, left_index=True, right_index=True, how="inner")
data = pd.merge(data, groups, left_index=True, right_index=True, how="inner")
data = data[data["group"] != "QC"]
groups_4 = data["group"].copy()
data["group"] = data["group"] == "HCC"

clinical_features = clinical.columns.tolist()
glycan_features = abundance.columns.tolist()

train_data, test_data = train_test_split(data, test_size=128, random_state=2, stratify=groups_4, shuffle=True)
feature_types = {"clinical": clinical_features, "glycan": glycan_features}

train_data.to_csv(snakemake.output["train_data"], index=True)
test_data.to_csv(snakemake.output["test_data"], index=True)
with open(snakemake.output["feature_types"], "w") as f:
    json.dump(feature_types, f)
