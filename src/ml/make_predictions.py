import json

import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

from hcc_fusion_model import HCCFusionClassifier


train_data = pd.read_csv(snakemake.input["train_data"], index_col=0)
test_data = pd.read_csv(snakemake.input["test_data"], index_col=0)
feature_types = json.load(open(snakemake.input["feature_types"]))

X_train = train_data.drop(columns=["group"])
y_train = train_data["group"]
X_test = test_data.drop(columns=["group"])
y_test = test_data["group"]

complex_model = HCCFusionClassifier(
    clinical_features=feature_types["clinical"], 
    glycan_features=feature_types["glycan"],
    random_state=42,
)
complex_model = CalibratedClassifierCV(complex_model, method="isotonic", cv=5)
complex_model.fit(X_train, y_train)

selected_features = ["H5N4F1", "H6N5F1S3", "H3N4F1"]
simple_model = make_pipeline(StandardScaler(), LogisticRegression(random_state=42))
simple_model.fit(X_train[selected_features], y_train)

y_pred_complex = complex_model.predict(X_test)
y_proba_complex = complex_model.predict_proba(X_test)[:, 1]
y_pred_simple = simple_model.predict(X_test[selected_features])
y_proba_simple = simple_model.predict_proba(X_test[selected_features])[:, 1]

prediction_df_complex = pd.DataFrame(
    {
        "model": "HCC Fusion",
        "target": y_test,
        "prediction": y_pred_complex,
        "probability": y_proba_complex,
    },
    index=X_test.index,
)
prediction_df_simple = pd.DataFrame(
    {
        "model": "HCC Slim",
        "target": y_test,
        "prediction": y_pred_simple,
        "probability": y_proba_simple,
    },
    index=X_test.index,
)
prediction_df = pd.concat([prediction_df_complex, prediction_df_simple])
prediction_df.to_csv(snakemake.output[0], index=True)
