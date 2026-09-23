import os
import sys
import numpy as np
import warnings

warnings.filterwarnings('ignore')



from config import LABEL_COLS, OUT_DIR
from data_utils import load_data, prepare_features_and_labels_per_label, split_data
from pipeline import run_experiment


def main():

    df = load_data()

    all_summary = []

    for label_name in LABEL_COLS:
        print(f"  Label: {label_name}")

        X, y, df_valid = prepare_features_and_labels_per_label(df, label_name)

        print(f"  Positives: {np.sum(y == 1)}/{len(y)} "
              f"({np.mean(y == 1)*100:.2f}%)")

        X_train, y_train, X_val, y_val, X_test, y_test = split_data(X, y, df_valid=df_valid)
        print(f"  Split: Train={len(X_train)} | Val={len(X_val)} | "
              f"Test={len(X_test)}")

        groups_to_keep = ['morgan', 'desc']
            
        exp_res = run_experiment(
            label=label_name,
            X_train=X_train, y_train=y_train,
            X_val=X_val, y_val=y_val,
            X_test=X_test, y_test=y_test,
            groups_to_keep=groups_to_keep,
            output_dir=OUT_DIR
        )
        all_summary.append(exp_res)

    print("Summary")

    for res in all_summary:
        lbl = res['label']
        best_m = res['best_model_name']
        auc = res['test_metrics']['ROC_AUC']
        f1 = res['test_metrics']['F1']
        pr_auc = res['test_metrics']['PR_AUC']
        print(f"  [{lbl:15s}] Best: {best_m:18s} | "
              f"Test ROC_AUC={auc:.4f} | PR_AUC={pr_auc:.4f} | F1={f1:.4f}")

    print(f"\n models & results saved to: {OUT_DIR}/")


if __name__ == "__main__":
    main()
