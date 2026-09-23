# =============================================================================
# Ghép metadata gốc (Pattern, Source, PubChem_ID, norm) trở lại SMILES đã curate.
# Inner join theo Drug_Name -> chỉ giữ drug survive MEHC.
# =============================================================================
import pandas as pd

INTERMEDIATE_DIR = "../../../../data/01_intermediate/phase1_curation/adrecs"
RAW = f"{INTERMEDIATE_DIR}/adrecs_chembl_smiles.csv"
CUR = f"{INTERMEDIATE_DIR}/mehc/refinement/post_refined_smiles.csv"
OUT = "../../../../data/02_curated_raw/adrecs_pneumotox_curated.csv"

def main():
    raw = pd.read_csv(RAW, encoding="utf-8-sig")
    cur = pd.read_csv(CUR, encoding="utf-8-sig")

    # 1 dòng / drug, giữ bản ghi đầu tiên
    raw_u = raw.drop_duplicates(subset=["Drug_Name"], keep="first")
    cur_u = cur.drop_duplicates(subset=["Drug_Name"], keep="first")

    # Inner join: dùng suffixes để tránh lỗi trùng tên cột SMILES
    merged = raw_u.merge(
        cur_u[["Drug_Name", "smiles"]].rename(columns={"smiles": "SMILES"}),
        on="Drug_Name", how="inner", suffixes=("_raw", "_cur")
    )

    # Lấy SMILES chuẩn từ MEHC (SMILES_cur), xoá SMILES_raw cũ
    if "SMILES_raw" in merged.columns:
        merged["SMILES"] = merged["SMILES_cur"]
        merged = merged.drop(columns=["SMILES_raw", "SMILES_cur"])

    # Đảm bảo cột SMILES nằm cuối (nếu raw chưa có)
    cols = list(raw.columns)
    if "SMILES" not in cols:
        cols.append("SMILES")
    merged = merged[cols]

    merged.to_csv(OUT, index=False, encoding="utf-8-sig")
    print(f"{OUT}: {len(merged)} rows")

if __name__ == "__main__":
    main()