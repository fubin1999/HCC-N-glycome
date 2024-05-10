import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin, ClassifierMixin
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.multiclass import unique_labels
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


class HCCFusionClassifier(BaseEstimator, TransformerMixin, ClassifierMixin):
    """The final model for HCC diagnosis.

    This model is a fusion of two models: 
    1. A Random Forest model trained on the clinical data;
    2. A SVM model trained on the glycomics data.

    The two models are trained separately and their predicted probabilities
    are combined using a meta-classifier (Logistic Regression) to make the final prediction.
    """

    def __init__(self, clinical_features, glycan_features, random_state=None):
        super().__init__()
        self.clinical_features = clinical_features
        self.glycan_features = glycan_features
        self.random_state = random_state

        self.clinical_model = RandomForestClassifier(n_estimators=1000, random_state=random_state)
        self.glycan_model = make_pipeline(StandardScaler(), SVC(C=5, gamma=0.01, probability=True, random_state=random_state))
        self.meta_classifier = LogisticRegression(random_state=random_state)

    def fit(self, X, y):
        self.classes_ = unique_labels(y)
        self.clinical_model.fit(X[self.clinical_features], y)
        self.glycan_model.fit(X[self.glycan_features], y)
        X_meta = self._get_X_meta(X)
        self.meta_classifier.fit(X_meta, y)
        return self

    def predict(self, X):
        check_is_fitted(self)
        X_meta = self._get_X_meta(X)
        return self.meta_classifier.predict(X_meta)

    def predict_proba(self, X):
        check_is_fitted(self)
        X_meta = self._get_X_meta(X)
        return self.meta_classifier.predict_proba(X_meta)

    def _get_X_meta(self, X):
        proba_clinical = self.clinical_model.predict_proba(X[self.clinical_features])[:, 1]
        proba_glycan = self.glycan_model.predict_proba(X[self.glycan_features])[:, 1]
        X_meta = np.vstack([proba_clinical, proba_glycan]).T
        return X_meta
