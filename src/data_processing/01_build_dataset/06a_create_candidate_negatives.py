# =============================================================================
# Lọc các row có TỔNG 4 nhãn = 0 -> candidate negatives (chờ chấm qua PubChem).
# =============================================================================
import pandas as pd

IN  = '../../../data/01_intermediate/phase2_assembly/05_merged_ctd_positive_pool.csv'
OUT = '../../../data/01_intermediate/phase2_assembly/06_candidate_negatives_raw.csv'
LABEL_COLS = ['Label_ILD', 'Label_Edema', 'Label_Pleural', 'Label_Embolism']

df = pd.read_csv(IN, encoding='utf-8')
neg = df[df[LABEL_COLS].sum(axis=1) == 0].copy()
neg.to_csv(OUT, index=False, encoding='utf-8')
print(f"{OUT}: {len(neg)} rows")