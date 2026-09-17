"""
12 model configs x 4 targets x 10 seeds.
Feature groups:
  - desc:           RDKit 2D descriptors only
  - morgan:         Morgan ECFP4 fingerprint only
  - desc+morgan:    Both combined

Classifiers: KNN, SVM, RandomForest, XGBoost
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from collections import defaultdict
from sklearn.model_selection import StratifiedKFold, GridSearchCV
import warnings
from rdkit import RDLogger

warnings.filterwarnings('ignore')
RDLogger.DisableLog('rdApp.*')



from config import LABEL_COLS, CV_FOLDS, get_feature_groups
from data_utils import load_data, filter_clean_negatives, split_data
from featurizer import featurize
from models import get_search_spaces
from evaluator import evaluate


SEEDS = list(range(10))         
N_JOBS = -3
FIXED_CV_SEED = 42               

FEATURE_CONFIGS = {
    "Desc":         ["desc"],
    "Morgan":       ["morgan"],
    "Desc+Morgan":  ["morgan", "desc"],
}
CLASSIFIERS = ["KNN", "SVM", "RandomForest", "XGBoost"]

CSV_PATH = "balanced_pulmonary_adr_dataset.csv"


def run_single(label_name, feature_groups, classifier_name, X, y, df_valid, seed, split_type="scaffold"):
    if split_type == "scaffold":
        X_train, y_train, X_val, y_val, X_test, y_test = split_data(
            X, y, df_valid=df_valid, random_state=seed
        )
    else:
        X_train, y_train, X_val, y_val, X_test, y_test = split_data(
            X, y, df_valid=None, random_state=seed
        )
    searches = get_search_spaces(
        groups_to_keep=feature_groups,
        y_train=y_train,
        random_state=FIXED_CV_SEED
    )

    pipe, param_grid = searches[classifier_name]
    cv = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=FIXED_CV_SEED)

    gs = GridSearchCV(
        pipe, param_grid, cv=cv, scoring='roc_auc',
        n_jobs=N_JOBS, refit=True, verbose=0
    )
    gs.fit(X_train, y_train)
    best_est = gs.best_estimator_

    # Predict on test set
    if hasattr(best_est, 'decision_function'):
        test_scores = best_est.decision_function(X_test)
    else:
        test_scores = best_est.predict_proba(X_test)[:, 1]
    test_pred = best_est.predict(X_test)

    metrics = evaluate(y_test, test_scores, test_pred)
    return metrics, (y_test, test_scores, test_pred)


def main():
    t0 = time.time()

    df = load_data(CSV_PATH)

    all_results = defaultdict(lambda: defaultdict(list))
    predictions_list = []

    for label_name in LABEL_COLS:
        print(f"  TARGET: {label_name}")
        df_clean = filter_clean_negatives(df, label_name)
        smiles_list = df_clean["SMILES"].tolist()
        X_full, valid_idx = featurize(smiles_list)
        X_full = np.nan_to_num(X_full.astype(np.float32))
        df_valid = df_clean.iloc[valid_idx].reset_index(drop=True)
        y = df_valid[label_name].values.astype(int)

        print(f"  Samples: {len(y)} | Pos: {y.sum()} | Neg: {len(y)-y.sum()}")
        print(f"  Feature matrix: {X_full.shape}")

        for split_type in ["scaffold", "random"]:
            for feat_name, feat_groups in FEATURE_CONFIGS.items():
                for clf_name in CLASSIFIERS:
                    config_name = f"{split_type} + {feat_name} + {clf_name}"
                    print(f"\n  >> {config_name}")

                    for seed in SEEDS:
                        metrics, preds_tuple = run_single(
                            label_name, feat_groups, clf_name,
                            X_full, y, df_valid, seed, split_type
                        )
                        all_results[label_name][config_name].append(metrics)
                        
                        y_true_seed, y_prob_seed, y_pred_seed = preds_tuple
                        df_pred = pd.DataFrame({
                            "Target": label_name,
                            "Split": split_type,
                            "Features": feat_name,
                            "Classifier": clf_name,
                            "Seed": seed,
                            "y_true": y_true_seed,
                            "y_prob": y_prob_seed,
                            "y_pred": y_pred_seed
                        })
                    predictions_list.append(df_pred)
                    print(f"    [Seed {seed}] ROC_AUC: {metrics['ROC_AUC']:.3f} | PR_AUC: {metrics['PR_AUC']:.3f} | F1: {metrics['F1']:.3f}")

                rocs = [m['ROC_AUC'] for m in all_results[label_name][config_name]]
                print(f"  -> Mean ROC_AUC: {np.mean(rocs):.3f} +/- {np.std(rocs):.3f}")

    print("  SUMMARY: MEAN +/- STD ACROSS 10 SEEDS")
    rows = []
    for label_name in LABEL_COLS:
        target = label_name.replace("Label_", "")
        for split_type in ["scaffold", "random"]:
            for feat_name in FEATURE_CONFIGS:
                for clf_name in CLASSIFIERS:
                    config_name = f"{split_type} + {feat_name} + {clf_name}"
                    metrics_list = all_results[label_name][config_name]
                    row = {
                        "Target": target,
                        "Split": split_type,
                        "Features": feat_name,
                        "Classifier": clf_name,
                    }
                    for metric_key in ["ROC_AUC", "PR_AUC", "F1", "Precision", "Recall", "Balanced_Acc"]:
                        vals = [m[metric_key] for m in metrics_list]
                        row[metric_key] = f"{np.mean(vals):.3f} +/- {np.std(vals):.3f}"
                    rows.append(row)

    df_summary = pd.DataFrame(rows)

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 200)
    print(df_summary.to_string(index=False))

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, "baseline_ablation_results.csv")
    df_summary.to_csv(out_csv, index=False)

    df_preds_all = pd.concat(predictions_list, ignore_index=True)
    out_preds_csv = os.path.join(out_dir, "baseline_ablation_predictions.csv")
    df_preds_all.to_csv(out_preds_csv, index=False)

    elapsed = time.time() - t0
    print(f"Total time: {elapsed/60:.1f} min")


if __name__ == "__main__":
    main()
