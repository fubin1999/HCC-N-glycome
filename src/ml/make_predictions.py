import json

import pandas as pd
import numpy as np

from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
)

from hcc_fusion_model import HCCFusionClassifier


train_data = pd.read_csv(snakemake.input["train_data"], index_col=0)
test_data = pd.read_csv(snakemake.input["test_data"], index_col=0)
feature_types = json.load(open(snakemake.input["feature_types"]))

X_train = train_data.drop(columns=["group"])
y_train = train_data["group"]
X_test = test_data.drop(columns=["group"])
y_test = test_data["group"]

model = HCCFusionClassifier(
    clinical_features=feature_types["clinical"], 
    glycan_features=feature_types["glycan"],
    random_state=42,
)
model.fit(X_train, y_train)

y_pred = model.predict(X_test)
y_proba = model.predict_proba(X_test)[:, 1]

metrics = {
    "accuracy": accuracy_score(y_test, y_pred),
    "precision": precision_score(y_test, y_pred),
    "recall": recall_score(y_test, y_pred),
    "f1": f1_score(y_test, y_pred),
    "roc_auc": roc_auc_score(y_test, y_proba),
}
json.dump(metrics, open(snakemake.output["metrics"], "w"))

prediction_df = pd.DataFrame(
    {
        "target": y_test,
        "prediction": y_pred,
        "probability": y_proba,
    },
    index=X_test.index,
)
prediction_df.to_csv(snakemake.output["predictions"], index=True)
