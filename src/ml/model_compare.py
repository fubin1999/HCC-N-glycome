from collections import namedtuple
import warnings

import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression, Perceptron, PassiveAggressiveClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.svm import SVC, LinearSVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF, DotProduct, Matern, WhiteKernel, RationalQuadratic
from sklearn.naive_bayes import GaussianNB, MultinomialNB, ComplementNB, BernoulliNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.decomposition import PCA
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier, ExtraTreesClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.preprocessing import StandardScaler

from hcc_fusion_model import HCCFusionClassifier

warnings.filterwarnings("ignore", module="sklearn")


# abundance = pd.read_csv("results/data/prepared/processed_abundance.csv")
# clinical = pd.read_csv("results/data/prepared/clinical.csv")
# groups = pd.read_csv("results/data/prepared/groups.csv")
abundance = pd.read_csv(snakemake.input["abundance"])
clinical = pd.read_csv(snakemake.input["clinical"])
groups = pd.read_csv(snakemake.input["groups"])

abundance = abundance.pivot(index="sample", columns="glycan", values="value")
abundance = pd.DataFrame(np.log2(abundance.values), columns=abundance.columns, index=abundance.index)
clinical = clinical.drop(["sex", "age"], axis=1)
clinical = clinical.set_index("sample")
groups = groups.set_index("sample")
data = pd.merge(abundance, clinical, left_index=True, right_index=True, how="inner")
data = pd.merge(data, groups, left_index=True, right_index=True, how="inner")
data = data[data["group"] != "QC"]
groups_4 = data["group"].copy()
data["group"] = data["group"] == "C"

clinical_features = clinical.columns.tolist()
glycan_features = abundance.columns.tolist()

train_data, test_data = train_test_split(data, test_size=128, random_state=42, stratify=groups_4, shuffle=True)

X_train = train_data.drop("group", axis=1)
y_train = train_data["group"]
X_test = test_data.drop("group", axis=1)
y_test = test_data["group"]


models = {
    "Logistic Regression": LogisticRegression(penalty=None, max_iter=1000),
    "Lasso Logistic Regression (C=1)": LogisticRegression(penalty="l1", C=1, solver="saga", max_iter=1000),
    "Lasso Logistic Regression (C=10)": LogisticRegression(penalty="l1", C=10, solver="saga", max_iter=1000),
    "Lasso Logistic Regression (C=0.1)": LogisticRegression(penalty="l1", C=0.1, solver="saga", max_iter=1000),
    "Ridge Logistic Regression (C=1)": LogisticRegression(penalty="l2", C=1, max_iter=1000),
    "Ridge Logistic Regression (C=10)": LogisticRegression(penalty="l2", C=10, max_iter=1000),
    "Ridge Logistic Regression (C=0.1)": LogisticRegression(penalty="l2", C=0.1, max_iter=1000),
    "Elastic Net (C=1)": LogisticRegression(penalty="elasticnet", l1_ratio=0.5, C=1, solver="saga", max_iter=1000),
    "Elastic Net (C=10)": LogisticRegression(penalty="elasticnet", l1_ratio=0.5, C=10, solver="saga", max_iter=1000),
    "Elastic Net (C=0.1)": LogisticRegression(penalty="elasticnet", l1_ratio=0.5, C=0.1, solver="saga", max_iter=1000),
    "Perceptron": Perceptron(max_iter=1000),
    "Passive-Aggressive (C=1)": PassiveAggressiveClassifier(C=1, max_iter=1000),
    "Passive-Aggressive (C=10)": PassiveAggressiveClassifier(C=10, max_iter=1000),
    "Passive-Aggressive (C=0.1)": PassiveAggressiveClassifier(C=0.1, max_iter=1000),
    "Linear Discriminant Analysis": LinearDiscriminantAnalysis(),
    "Quadratic Discriminant Analysis": QuadraticDiscriminantAnalysis(),
    "Linear SVM (C=1)": LinearSVC(C=1, max_iter=1000),
    "Linear SVM (C=10)": LinearSVC(C=10, max_iter=1000),
    "Linear SVM (C=0.1)": LinearSVC(C=0.1, max_iter=1000),
    "RBF SVM (C=1, gamma=0.01)": SVC(C=1, gamma=0.01, probability=True, max_iter=1000),
    "RBF SVM (C=10, gamma=0.01)": SVC(C=10, gamma=0.01, probability=True, max_iter=1000),
    "RBF SVM (C=0.1, gamma=0.01)": SVC(C=0.1, gamma=0.01, probability=True, max_iter=1000),
    "RBF SVM (C=1, gamma=0.1)": SVC(C=1, gamma=0.1, probability=True, max_iter=1000),
    "RBF SVM (C=10, gamma=0.1)": SVC(C=10, gamma=0.1, probability=True, max_iter=1000),
    "RBF SVM (C=0.1, gamma=0.1)": SVC(C=0.1, gamma=0.1, probability=True, max_iter=1000),
    "Polynomial SVM (C=1, degree=2)": SVC(C=1, kernel="poly", degree=2, probability=True, max_iter=1000),
    "Polynomial SVM (C=10, degree=2)": SVC(C=10, kernel="poly", degree=2, probability=True, max_iter=1000),
    "Polynomial SVM (C=0.1, degree=2)": SVC(C=0.1, kernel="poly", degree=2, probability=True, max_iter=1000),
    "Polynomial SVM (C=1, degree=3)": SVC(C=1, kernel="poly", degree=3, probability=True, max_iter=1000),
    "Polynomial SVM (C=10, degree=3)": SVC(C=10, kernel="poly", degree=3, probability=True, max_iter=1000),
    "Polynomial SVM (C=0.1, degree=3)": SVC(C=0.1, kernel="poly", degree=3, probability=True, max_iter=1000),
    "Polynomial SVM (C=1, degree=4)": SVC(C=1, kernel="poly", degree=4, probability=True, max_iter=1000),
    "Polynomial SVM (C=10, degree=4)": SVC(C=10, kernel="poly", degree=4, probability=True, max_iter=1000),
    "Polynomial SVM (C=0.1, degree=4)": SVC(C=0.1, kernel="poly", degree=4, probability=True, max_iter=1000),
    "K-Nearest Neighbors (k=5)": KNeighborsClassifier(n_neighbors=5),
    "K-Nearest Neighbors (k=10)": KNeighborsClassifier(n_neighbors=10),
    "K-Nearest Neighbors (k=15)": KNeighborsClassifier(n_neighbors=20),
    "Gaussian Process (RBF kernel)": GaussianProcessClassifier(kernel=1.0 * RBF(1.0)),
    "Gaussian Process (Dot Product kernel)": GaussianProcessClassifier(kernel=1.0 * DotProduct(1.0)),
    "Gaussian Process (Matern kernel)": GaussianProcessClassifier(kernel=1.0 * Matern(1.0)),
    "Gaussian Process (RBF kernel + White Noise)": GaussianProcessClassifier(kernel=1.0 * RBF(1.0) + WhiteKernel(1.0)),
    "Gaussian Process (Rantional Quadratic kernel)": GaussianProcessClassifier(kernel=1.0 * RationalQuadratic(1.0)),
    "Gaussian Naive Bayes": GaussianNB(),
    "Multinomial Naive Bayes": MultinomialNB(),
    "Complement Naive Bayes": ComplementNB(),
    "Bernoulli Naive Bayes": BernoulliNB(),
    "PCA + Decision Tree (max_depth=3)": make_pipeline(PCA(), DecisionTreeClassifier(max_depth=3)),
    "PCA + Decision Tree (max_depth=5)": make_pipeline(PCA(), DecisionTreeClassifier(max_depth=5)),
    "PCA + Decision Tree (max_depth=7)": make_pipeline(PCA(), DecisionTreeClassifier(max_depth=7)),
    "PCA + Decision Tree (max_depth=9)": make_pipeline(PCA(), DecisionTreeClassifier(max_depth=9)),
    "Gradient Tree Boosting (max_depth=3)": GradientBoostingClassifier(max_depth=3),
    "Gradient Tree Boosting (max_depth=5)": GradientBoostingClassifier(max_depth=5),
    "Gradient Tree Boosting (max_depth=7)": GradientBoostingClassifier(max_depth=7),
    "Gradient Tree Boosting (max_depth=9)": GradientBoostingClassifier(max_depth=9),
    "Random Forest (max_depth=3)": RandomForestClassifier(max_depth=3),
    "Random Forest (max_depth=5)": RandomForestClassifier(max_depth=5),
    "Random Forest (max_depth=7)": RandomForestClassifier(max_depth=7),
    "Random Forest (max_depth=9)": RandomForestClassifier(max_depth=9),
    "Extra Trees (max_depth=3)": ExtraTreesClassifier(max_depth=3),
    "Extra Trees (max_depth=5)": ExtraTreesClassifier(max_depth=5),
    "Extra Trees (max_depth=7)": ExtraTreesClassifier(max_depth=7),
    "Extra Trees (max_depth=9)": ExtraTreesClassifier(max_depth=9),
    "AdaBoost (max_depth=3)": AdaBoostClassifier(DecisionTreeClassifier(max_depth=3)),
    "AdaBoost (max_depth=5)": AdaBoostClassifier(DecisionTreeClassifier(max_depth=5)),
    "AdaBoost (max_depth=7)": AdaBoostClassifier(DecisionTreeClassifier(max_depth=7)),
    "AdaBoost (max_depth=9)": AdaBoostClassifier(DecisionTreeClassifier(max_depth=9)),
}

no_scale_keywords = ["Naive Bayes", "Tree", "Forest", "Boost"]
for model_name, model in models.items():
    model.random_state = 42
    if any(keyword in model_name for keyword in no_scale_keywords):
        models[model_name] = model
    else:
        models[model_name] = make_pipeline(StandardScaler(), model)
models["HCC Fusion Classifier"] = HCCFusionClassifier(clinical_features, glycan_features, random_state=42)

Record = namedtuple("Record", ["acc_mean", "auc_mean", "f1_mean", "acc_std", "auc_std", "f1_std"])
results = {}
for model_name, model in models.items():
    print(f"===== {model_name} =====")
    accs = cross_val_score(model, X_train, y_train, cv=10, scoring="accuracy")
    acc_mean = accs.mean()
    acc_std = accs.std()
    print(f"Accuracy: {acc_mean:.3f} +- {acc_std:.3f}")
    aucs = cross_val_score(model, X_train, y_train, cv=10, scoring="roc_auc")
    auc_mean = aucs.mean()
    auc_std = aucs.std()
    print(f"ROC AUC: {auc_mean:.3f} +- {auc_std:.3f}")
    f1s = cross_val_score(model, X_train, y_train, cv=10, scoring="f1")
    f1_mean = f1s.mean()
    f1_std = f1s.std()
    print(f"F1-score: {f1_mean:.3f} +- {f1_std:.3f}")
    results[model_name] = Record(acc_mean, auc_mean, f1_mean, acc_std, auc_std, f1_std)
    print()

results = pd.DataFrame(results).T
results.columns = Record._fields
# results.to_csv("results/data/ml/model_comparison.csv")
results.to_csv(snakemake.output[0])
