import numpy as np
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from evaluator import evaluate, print_metrics_block
from models import get_search_spaces
from config import CV_FOLDS, OUT_DIR
import joblib
import os

def run_experiment(label, X_train, y_train, X_val, y_val, X_test, y_test,
                   groups_to_keep, output_dir=OUT_DIR, models_to_run=None, random_state=42):
    searches = get_search_spaces(groups_to_keep, y_train=y_train, random_state=random_state)
    
    if models_to_run:
        searches = {k: v for k, v in searches.items() if k in models_to_run}

    cv = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=random_state)

    best_val_auc = -1
    best_model_name = None
    best_pipe = None
    best_test_metrics = None
    best_val_metrics = None
    best_params = None
    results = {}

    for model_name, (pipe, param_grid) in searches.items():
        gs = GridSearchCV(
            pipe, param_grid, cv=cv, scoring='roc_auc',
            n_jobs=4, refit=True, verbose=0
        )
        gs.fit(X_train, y_train)
        best_est = gs.best_estimator_

        # ---- Validation ----
        if hasattr(best_est, 'decision_function'):
            val_scores = best_est.decision_function(X_val)
        else:
            val_scores = best_est.predict_proba(X_val)[:, 1]
        val_pred = best_est.predict(X_val)
        val_metrics = evaluate(y_val, val_scores, val_pred)

        # ---- Test ----
        if hasattr(best_est, 'decision_function'):
            test_scores = best_est.decision_function(X_test)
        else:
            test_scores = best_est.predict_proba(X_test)[:, 1]
        test_pred = best_est.predict(X_test)
        test_metrics = evaluate(y_test, test_scores, test_pred)

        results[model_name] = {
            'best_estimator': best_est,
            'val_metrics': val_metrics,
            'test_metrics': test_metrics,
            'best_params': gs.best_params_,
            'cv_mean': gs.cv_results_['mean_test_score'][gs.best_index_],
            'cv_std': gs.cv_results_['std_test_score'][gs.best_index_],
        }

        print(f"\n  ── {model_name} ──  best_params={gs.best_params_}")
        print(f"    CV (train)   ROC_AUC={gs.cv_results_['mean_test_score'][gs.best_index_]:.3f}±"
              f"{gs.cv_results_['std_test_score'][gs.best_index_]:.3f}")
        print_metrics_block("Validation", val_metrics)
        print_metrics_block("Test", test_metrics)

        if val_metrics['ROC_AUC'] > best_val_auc:
            best_val_auc = val_metrics['ROC_AUC']
            best_model_name = model_name
            best_pipe = best_est
            best_test_metrics = test_metrics
            best_val_metrics = val_metrics
            best_params = gs.best_params_


    champion_name = best_model_name
    champion_val_auc = best_val_auc
    champion_test_metrics = best_test_metrics

    # Save best model
    joblib.dump(best_pipe,
                os.path.join(output_dir, f"best_model_{label}.joblib"))

    print(f"\n  >>> Best model for '{label}': {champion_name} "
          f"(Val-ROC-AUC={champion_val_auc:.3f})")
    print_metrics_block("Test-final", champion_test_metrics)

    return {
        'label': label,
        'groups': groups_to_keep,
        'best_model_name': champion_name,
        'test_metrics': champion_test_metrics,
        'all_results': results,
    }