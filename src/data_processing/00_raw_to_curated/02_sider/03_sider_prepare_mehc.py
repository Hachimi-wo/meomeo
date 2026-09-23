# =============================================================================
# 03_sider_prepare_mehc.py
# Format input cho MEHC: SMILES + index (1-based).
# 'index' là join key để merge metadata lại ở bước 05.
# =============================================================================
import pandas as pd

INTERMEDIATE_DIR = "../../../../data/01_intermediate/phase1_curation/sider"
IN  = f'{INTERMEDIATE_DIR}/sider_filtered_nomehc.csv'
OUT = f'{INTERMEDIATE_DIR}/respiratory_labeled_for_curation.csv'

def main():
    df = pd.read_csv(IN)
    out = df[['SMILES']].copy()
    out.insert(1, 'index', range(1, len(out) + 1))
    out.to_csv(OUT, index=False, encoding='utf-8-sig')

if __name__ == "__main__":
    main()