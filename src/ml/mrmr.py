import json

import pandas as pd
import numpy as np

from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression
from mrmr import mrmr_classif


train_data = pd.read_csv(snakemake.input[0], index_col=0)
X_train = train_data.drop("group", axis=1)
y_train = train_data["group"]

feature_types = json.load(open(snakemake.input[1]))
glycan_features = feature_types["glycan"]
X_train_glycan = X_train[glycan_features]

score_means = np.empty(len(glycan_features))
score_stds = np.empty(len(glycan_features))

ranks = mrmr_classif(X_train_glycan, y_train, K=len(glycan_features), show_progress=False)

for i in range(len(glycan_features)):
    selected_glycans = ranks[:i+1]
    X_train_subset = X_train[selected_glycans + ["AFP"]]
    pipe = make_pipeline(StandardScaler(), LogisticRegression())
    scores = cross_val_score(pipe, X_train_subset, y_train, cv=10, scoring="roc_auc")
    score_means[i] = scores.mean()
    score_stds[i] = scores.std()

result_df = pd.DataFrame({
    "n_features": list(range(1, len(glycan_features)+1)),
    "glycan": ranks,
    "score_mean": score_means,
    "score_std": score_stds
})
result_df.to_csv(snakemake.output[0], index=False)
