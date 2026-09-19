import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, RobustScaler, MaxAbsScaler
from sklearn.compose import ColumnTransformer
from sklearn.feature_selection import VarianceThreshold
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.calibration import CalibratedClassifierCV
from xgboost import XGBClassifier
from imblearn.pipeline import Pipeline as ImbPipeline
from imblearn.under_sampling import RandomUnderSampler
from config import N_ESTIMATORS, get_feature_groups

def _compute_pos_weight(y_train):
    n_neg = np.sum(y_train == 0)
    n_pos = np.sum(y_train == 1)
    if n_pos == 0:
        return 1.0
    return float(n_neg) / float(n_pos)


def get_search_spaces(groups_to_keep, y_train=None, random_state=42):
    groups = get_feature_groups()
    transformers = []
    
    if 'morgan' in groups_to_keep and 'morgan' in groups:
        start, end = groups['morgan']
        cols = list(range(start, end))
        morgan_pipe = Pipeline([
            ('imputer', SimpleImputer(strategy='constant', fill_value=0.0)),
            ('varth', VarianceThreshold(threshold=0.01)),
            ('scaler', MaxAbsScaler())
        ])
        transformers.append(('morgan', morgan_pipe, cols))

    if 'desc' in groups_to_keep and 'desc' in groups:
        start, end = groups['desc']
        from featurizer import _init_descriptors, _DESC_NAMES
        _init_descriptors()
        discrete_cols = []
        continuous_cols = []
        for i, name in enumerate(_DESC_NAMES):
            idx = start + i
            if "Count" in name or "Num" in name or "fr_" in name:
                discrete_cols.append(idx)
            else:
                continuous_cols.append(idx)
        
        if discrete_cols:
            disc_pipe = Pipeline([
                ('imputer', SimpleImputer(strategy='median')),
                ('scaler', RobustScaler())
            ])
            transformers.append(('desc_disc', disc_pipe, discrete_cols))
            
        if continuous_cols:
            cont_pipe = Pipeline([
                ('imputer', SimpleImputer(strategy='mean')),
                ('scaler', RobustScaler())
            ])
            transformers.append(('desc_cont', cont_pipe, continuous_cols))

    if not transformers:
        raise ValueError("No groups selected.")
    
    col_transformer = ColumnTransformer(transformers)

    common_steps = [
        ('preprocessor', col_transformer),
    ]

    searches = {}
    spw = _compute_pos_weight(y_train) if y_train is not None else 1.0

    # KNN (RandomUnderSampler)
    pipe_knn = ImbPipeline([
        ('preprocessor', col_transformer),
        ('rus', RandomUnderSampler(random_state=random_state)),
        ('clf', CalibratedClassifierCV(
            KNeighborsClassifier(weights='distance'),
            method='sigmoid',
            ensemble=False
        ))
    ])
    param_knn = {
        'clf__estimator__n_neighbors': list(range(3, 15)),
        'clf__estimator__p': [1, 2],
    }
    searches['KNN'] = (pipe_knn, param_knn)

    # SVM
    pipe_svm = Pipeline(common_steps + [
        ('clf', CalibratedClassifierCV(
            SVC(kernel='rbf', random_state=random_state),
            method='sigmoid',
            ensemble=False
        ))
    ])
    param_svm = {
        'clf__estimator__C': [0.001, 0.01, 0.1, 1, 10, 100, 1000],
        'clf__estimator__gamma': [0.001, 0.01, 0.1, 1, 10, 100, 1000],
        'clf__estimator__class_weight': [None, 'balanced'],
    }
    searches['SVM'] = (pipe_svm, param_svm)

    # Random Forest
    pipe_rf = Pipeline(common_steps + [
        ('clf', CalibratedClassifierCV(
            RandomForestClassifier(
                n_estimators=N_ESTIMATORS,
                random_state=random_state,
                n_jobs=1,
                class_weight='balanced'
            ),
            method='isotonic',
            ensemble=False
        ))
    ])
    param_rf = {
        'clf__estimator__max_depth': [3, 5, 7, 9],
        'clf__estimator__min_samples_split': [2, 4, 6, 8, 10],
        'clf__estimator__max_features': [0.2, 0.4, 0.6, 0.8],
    }
    searches['RandomForest'] = (pipe_rf, param_rf)

    # XGBoost
    pipe_xgb = Pipeline(common_steps + [
        ('clf', CalibratedClassifierCV(
            XGBClassifier(
                n_estimators=N_ESTIMATORS,
                random_state=random_state,
                n_jobs=1,
                eval_metric='logloss',
                scale_pos_weight=spw
            ),
            method='isotonic',
            ensemble=False
        ))
    ])
    param_xgb = {
        'clf__estimator__learning_rate': [0.0001, 0.001, 0.01, 0.1, 0.2, 0.3],
        'clf__estimator__max_depth': [3, 5, 7, 9],
        'clf__estimator__colsample_bytree': [0.3, 0.5, 0.7, 0.9],
    }
    searches['XGBoost'] = (pipe_xgb, param_xgb)

    return searches
