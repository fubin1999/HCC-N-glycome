from collections import defaultdict, namedtuple
from itertools import product

import pandas as pd
import numpy as np

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score

from tqdm import tqdm

# Read data
train_data = pd.read_csv(snakemake.input[0])
X_train = train_data.drop(columns=["group"])
y_train = train_data["group"]

# Read cluster data
cluster_df = pd.read_csv(snakemake.input[1])
clusters = defaultdict(list)
for _, row in cluster_df.iterrows():
    clusters[row["cluster"]].append(row["glycan"])

# Pick one feature from each cluster, all combinations
feature_combinations = []
for features in product(*clusters.values()):
    feature_combinations.append(list(features))

# Train a model for each feature combination
Record = namedtuple("Record", "features accuracy f1 roc_auc")
results: list[Record] = []
for features in tqdm(feature_combinations):
    features_str = "+".join(features)
    X_train_subset = X_train[features]
    model = make_pipeline(StandardScaler(), LogisticRegression(random_state=42))
    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    for train_index, test_index in skf.split(X_train_subset, y_train):
        model.fit(X_train_subset.iloc[train_index], y_train.iloc[train_index])
        y_pred = model.predict(X_train_subset.iloc[test_index])
        y_proba = model.predict_proba(X_train_subset.iloc[test_index])[:, 1]
        accuracy = accuracy_score(y_train.iloc[test_index], y_pred)
        f1 = f1_score(y_train.iloc[test_index], y_pred)
        roc_auc = roc_auc_score(y_train.iloc[test_index], y_proba)
        results.append(Record(features_str, accuracy, f1, roc_auc))

# Save results
results_df = pd.DataFrame(results, columns=Record._fields)
results_df.to_csv(snakemake.output[0], index=False)
