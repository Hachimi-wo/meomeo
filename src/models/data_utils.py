import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from config import (
    CSV_PATH, SMILES_COL, LABEL_COLS,
    FRAC_TRAIN, FRAC_VAL, FRAC_TEST
)
from featurizer import featurize


def load_data(csv_file=None):
    if csv_file is None:
        csv_file = CSV_PATH

    df = pd.read_csv(csv_file)
    
    n_before = len(df)
    df = df.drop_duplicates(subset=[SMILES_COL], keep='first').reset_index(drop=True)
    n_after = len(df)
    if n_before > n_after:
        print(f"Dropped {n_before - n_after} duplicate SMILES.")
    
    return df


def filter_clean_negatives(df, target_label):
    """
      - Positive: samples with target_label = 1
      - Clean Negative: samples with ALL 4 labels = 0
      - Excluded: negative samples in target_label but positive in other labels (ambiguous)

    Returns:
        df_clean: Filtered DataFrame (only pos + clean neg)
    """
    pos_mask = df[target_label] == 1
    all_neg_mask = df[LABEL_COLS].sum(axis=1) == 0

    keep_mask = pos_mask | all_neg_mask
    df_clean = df[keep_mask].reset_index(drop=True)

    n_pos = pos_mask.sum()
    n_clean_neg = all_neg_mask.sum()
    n_excluded = len(df) - keep_mask.sum()

    print(f"  Clean Negative Filter [{target_label}]:")
    print(f"    Pos: {n_pos} | Clean Neg: {n_clean_neg}  Excluded: {n_excluded}")
    print(f"    Size: {len(df)} -> {len(df_clean)}")

    return df_clean


def prepare_features_and_labels_per_label(df, target_label):
    df_clean = filter_clean_negatives(df, target_label)

    smiles_list = df_clean[SMILES_COL].tolist()
    X, valid_idx = featurize(smiles_list)
    df_valid = df_clean.iloc[valid_idx].reset_index(drop=True)

    X = np.asarray(X, dtype=np.float32)
    X[np.isinf(X)] = 0.0

    y = df_valid[target_label].values.astype(int)

    return X, y, df_valid


def prepare_features_and_labels(df):
    smiles_list = df[SMILES_COL].tolist()
    print(f"Featurizing {len(smiles_list)} SMILES...")
    X, valid_idx = featurize(smiles_list)
    X = np.asarray(X, dtype=np.float32)
    X[np.isinf(X)] = 0.0
    print(f"Valid SMILES: {len(valid_idx)}/{len(smiles_list)}")
    print(f"Feature matrix: {X.shape}")

    df_valid = df.iloc[valid_idx].reset_index(drop=True)
    labels_dict = {}
    for col in LABEL_COLS:
        if col in df_valid.columns:
            labels_dict[col] = df_valid[col].values.astype(int)

    return X, labels_dict, df_valid


def split_data(X, y, df_valid=None, random_state=42):
    if df_valid is None:
        val_test_size = FRAC_VAL + FRAC_TEST
        test_ratio_of_val_test = FRAC_TEST / val_test_size

        X_train, X_temp, y_train, y_temp = train_test_split(
            X, y, test_size=val_test_size, random_state=random_state, stratify=y
        )
        X_val, X_test, y_val, y_test = train_test_split(
            X_temp, y_temp, test_size=test_ratio_of_val_test,
            random_state=random_state, stratify=y_temp
        )
        return X_train, y_train, X_val, y_val, X_test, y_test

    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from collections import defaultdict
    import numpy as np

    #Scaffold-base split
    smiles_list = df_valid[SMILES_COL].tolist()
    scaffolds = defaultdict(list)
    
    for i, sm in enumerate(smiles_list):
        try:
            mol = Chem.MolFromSmiles(sm)
            if mol:
                scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
                scaffolds[scaffold].append(i)
            else:
                scaffolds[''].append(i)
        except:
            scaffolds[''].append(i)
            
    sorted_scaffolds = sorted(scaffolds.keys(), key=lambda k: (len(scaffolds[k]), k), reverse=True)
    scaffold_sets = [scaffolds[k] for k in sorted_scaffolds]
    
    rng = np.random.RandomState(random_state)
    rng.shuffle(scaffold_sets)
    
    train_idx, val_idx, test_idx = [], [], []
    train_cutoff = int(FRAC_TRAIN * len(smiles_list))
    val_cutoff = train_cutoff + int(FRAC_VAL * len(smiles_list))
    
    for scaffold_set in scaffold_sets:
        if len(train_idx) + len(scaffold_set) <= train_cutoff:
            train_idx.extend(scaffold_set)
        elif len(train_idx) + len(val_idx) + len(scaffold_set) <= val_cutoff:
            val_idx.extend(scaffold_set)
        else:
            test_idx.extend(scaffold_set)
            
    return X[train_idx], y[train_idx], X[val_idx], y[val_idx], X[test_idx], y[test_idx]
