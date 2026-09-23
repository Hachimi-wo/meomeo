# =============================================================================
# Ghép metadata (STITCH_ID, Drug_Name, Patterns) trở lại file curated qua `index`.
# Ghi đè sider_respiratory_adrs.csv bằng bản cuối cùng.
# =============================================================================
import pandas as pd

INTERMEDIATE_DIR = "../../../../data/01_intermediate/phase1_curation/sider"
MEHC_FILE = f"{INTERMEDIATE_DIR}/mehc/refinement/post_refined_smiles.csv"
RAW_FILE  = f"{INTERMEDIATE_DIR}/sider_filtered_nomehc_labeled.csv"
OUT_FILE  = "../../../../data/02_curated_raw/sider_respiratory_adrs.csv"

def main():
    df_mehc = pd.read_csv(MEHC_FILE)
    if "smiles" in df_mehc.columns:
        df_mehc = df_mehc.rename(columns={"smiles": "SMILES"})

    df_raw = pd.read_csv(RAW_FILE)
    df_raw['index'] = range(1, len(df_raw) + 1)

    # Left join giữ nguyên các dòng curated, chỉ lấy metadata tương ứng
    df_merged = pd.merge(
        df_mehc,
        df_raw[['index', 'STITCH_ID', 'Drug_Name', 'Respiratory_Patterns']],
        on='index', how='left',
    )

    cols = (['STITCH_ID', 'Drug_Name', 'SMILES', 'Respiratory_Patterns'])
    df_final = df_merged[cols]
    df_final.to_csv(OUT_FILE, index=False, encoding='utf-8-sig')
    print(f"{OUT_FILE}: {len(df_final)} rows")

if __name__ == "__main__":
    main()