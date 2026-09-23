from sklearn.metrics import (
    roc_auc_score, f1_score, precision_score,
    recall_score, balanced_accuracy_score,
    average_precision_score
)

def evaluate(y_true, y_score, y_pred=None):
    if y_pred is None:
        y_pred = (y_score >= 0.5).astype(int)
    return {
        "ROC_AUC": roc_auc_score(y_true, y_score),
        "PR_AUC": average_precision_score(y_true, y_score),
        "F1": f1_score(y_true, y_pred, zero_division=0),
        "Precision": precision_score(y_true, y_pred, zero_division=0),
        "Recall": recall_score(y_true, y_pred, zero_division=0),
        "Balanced_Acc": balanced_accuracy_score(y_true, y_pred),
    }

def print_metrics_block(title, metrics):
    print(f"    {title:12s} | " + " | ".join(f"{k}={v:.3f}" for k, v in metrics.items()))