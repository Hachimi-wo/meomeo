# =============================================================================
# 04_sider_run_mehc.py
# Chạy MEHC. Copy output -> sider_respiratory_adrs.csv
# =============================================================================
import os, shutil
import pandas as pd
from mehc_curation.refinement import RefinementStage

INTERMEDIATE_DIR = "../../../../data/01_intermediate/phase1_curation/sider"
INPUT = f"{INTERMEDIATE_DIR}/respiratory_labeled_for_curation.csv"
OUT_DIR = f"{INTERMEDIATE_DIR}/mehc"
FINAL = f"{INTERMEDIATE_DIR}/sider_respiratory_adrs.csv"

def main():
    df = pd.read_csv(INPUT)

    if "SMILES" in df.columns:
        df.rename(columns={"SMILES": "smiles"}, inplace=True)

    cols = ['smiles'] + [col for col in df.columns if col != 'smiles']
    df = df[cols]

    df = df.dropna(subset=['smiles']).reset_index(drop=True)

    try:
        RefinementStage(df).complete_refinement(
            output_dir=OUT_DIR, get_report=True, n_cpu=1,
        )
        shutil.copy(os.path.join(OUT_DIR, "refinement", "post_refined_smiles.csv"),
                    FINAL)
        print(f"Curated -> {FINAL}")
    except Exception as e:
        print(f"[Error] {e}")

if __name__ == "__main__":
    main()