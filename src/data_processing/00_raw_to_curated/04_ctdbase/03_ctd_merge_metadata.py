# =============================================================================
# 03_ctd_merge_metadata.py
# -----------------------------------------------------------------------------
# Description:
#   Ghép cột Pattern từ need2find.csv trở lại với file smiles đã qua MEHC.
# =============================================================================

import pandas as pd

INTERMEDIATE_DIR = "../../../../data/01_intermediate/phase1_curation/ctdbase"
RAW = "../../../../data/00_raw_data/ctdbase/need2find.csv"
CUR = f"{INTERMEDIATE_DIR}/mehc/refinement/post_refined_smiles.csv"
OUT = "../../../../data/02_curated_raw/ctd_mechanism_markers.csv"

def main():
    raw = pd.read_csv(RAW, encoding="utf-8-sig")
    cur = pd.read_csv(CUR, encoding="utf-8-sig")

    if "SMILES" not in cur.columns:
        if "smiles" in cur.columns:
            cur.rename(columns={"smiles": "SMILES"}, inplace=True)
            
    # Inner join on Drug_Name
    merged = cur.merge(raw[["Drug_Name", "Pattern"]], on="Drug_Name", how="inner")
    
    # Ensure correct column order: SMILES, Drug_Name, Pattern
    merged = merged[["SMILES", "Drug_Name", "Pattern"]]
    
    # Drop duplicates just in case
    merged = merged.drop_duplicates()

    merged.to_csv(OUT, index=False, encoding="utf-8-sig")
    print(f"Created file {OUT} with {len(merged)} rows (Pattern restored).")

if __name__ == "__main__":
    main()
