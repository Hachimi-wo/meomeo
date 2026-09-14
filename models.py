import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, RobustScaler, MaxAbsScaler
from sklearn.compose import ColumnTransformer
from sklearn.feature_selection import VarianceThreshold
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
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
    for g in groups_to_keep:
        if g not in groups:
            continue
        start, end = groups[g]
        cols = list(range(start, end))
        if g == 'desc':
            scaler = RobustScaler()
        else:
            scaler = Pipeline([
                ('varth', VarianceThreshold(threshold=0.01)),
                ('scaler', MaxAbsScaler())
            ])
        transformers.append((f'scaler_{g}', scaler, cols))

    if not transformers:
        raise ValueError("No groups selected.")
    col_transformer = ColumnTransformer(transformers)

    common_steps = [
        ('imputer', SimpleImputer(strategy='mean')),
        ('scaler', col_transformer),
    ]

    searches = {}
    spw = _compute_pos_weight(y_train) if y_train is not None else 1.0

    # KNN (RandomUnderSampler)
    pipe_knn = ImbPipeline([
        ('imputer', SimpleImputer(strategy='mean')),
        ('rus', RandomUnderSampler(random_state=random_state)),
        ('scaler', col_transformer),
        ('clf', KNeighborsClassifier(weights='distance'))
    ])
    param_knn = {
        'clf__n_neighbors': list(range(3, 15)),
        'clf__p': [1, 2],
    }
    searches['KNN'] = (pipe_knn, param_knn)

    # SVM
    pipe_svm = Pipeline(common_steps + [
        ('clf', SVC(kernel='rbf', random_state=random_state, probability=False))
    ])
    param_svm = {
        'clf__C': [0.001, 0.01, 0.1, 1, 10, 100, 1000],
        'clf__gamma': [0.001, 0.01, 0.1, 1, 10, 100, 1000],
        'clf__class_weight': [None, 'balanced'],
    }
    searches['SVM'] = (pipe_svm, param_svm)

    # Random Forest
    pipe_rf = Pipeline(common_steps + [
        ('clf', RandomForestClassifier(
            n_estimators=N_ESTIMATORS,
            random_state=random_state,
            n_jobs=1,
            class_weight='balanced'
        ))
    ])
    param_rf = {
        'clf__max_depth': [3, 5, 7, 9],
        'clf__min_samples_split': [2, 4, 6, 8, 10],
        'clf__max_features': [0.2, 0.4, 0.6, 0.8],
    }
    searches['RandomForest'] = (pipe_rf, param_rf)

    # XGBoost
    pipe_xgb = Pipeline(common_steps + [
        ('clf', XGBClassifier(
            n_estimators=N_ESTIMATORS,
            random_state=random_state,
            n_jobs=1,
            eval_metric='logloss',
            scale_pos_weight=spw
        ))
    ])
    param_xgb = {
        'clf__learning_rate': [0.0001, 0.001, 0.01, 0.1, 0.2, 0.3],
        'clf__max_depth': [3, 5, 7, 9],
        'clf__colsample_bytree': [0.3, 0.5, 0.7, 0.9],
    }
    searches['XGBoost'] = (pipe_xgb, param_xgb)

    return searches
