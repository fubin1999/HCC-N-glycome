from pathlib import Path

import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import matplotlib.pyplot as plt
import shap


train_data = pd.read_csv(snakemake.input[0], index_col=0)
test_data = pd.read_csv(snakemake.input[1], index_col=0)

features = ["H5N4F1", "H3N4F1", "H6N5F1S3", "H4N5F1S1", "H7N6F1S3", "H5N4F1S1"]

X_train = train_data[features]
y_train = train_data["group"]
X_test = test_data[features]
y_test = test_data["group"]

simple_model = make_pipeline(StandardScaler(), LogisticRegression(random_state=42))
simple_model.fit(X_train, y_train)

y_proba = simple_model.predict_proba(X_test)[:, 1]
indice = np.where((y_proba > 0.95) | (y_proba < 0.05) | ((y_proba > 0.45) & (y_proba < 0.55)))[0]

explainer = shap.Explainer(lambda x: simple_model.predict_proba(x)[:, 1], X_train)
shap_values = explainer(X_test)

shap.summary_plot(shap_values, show=False)
plt.savefig(snakemake.output[0], dpi=300)

waterfall_dirpath = Path(snakemake.output[1])
waterfall_dirpath.mkdir(parents=True, exist_ok=True)
for i in indice:
    sample_name = y_test.index[i]
    file_path = str(waterfall_dirpath / f"{sample_name}.pdf")
    fig = plt.figure()
    shap.waterfall_plot(shap_values[i], show=False)
    plt.gcf().set_size_inches(10, 6)
    plt.tight_layout()
    plt.savefig(file_path, dpi=300)
