# =============================================================================
# Chạy MEHC 
# =============================================================================

import os
import pandas as pd
from mehc_curation.refinement import RefinementStage

INTERMEDIATE_DIR = "../../../../data/01_intermediate/phase1_curation/drugcentral"
INPUT = f"{INTERMEDIATE_DIR}/drugcentralonly.csv"
OUT_DIR = f"{INTERMEDIATE_DIR}/mehc_output"
FINAL = "../../../../data/02_curated_raw/drugcentral_smiles.csv"

def main():
    df = pd.read_csv(INPUT)
    if "SMILES" in df.columns:
        df.rename(columns={"SMILES": "smiles"}, inplace=True)
    cols = ['smiles'] + [col for col in df.columns if col != 'smiles']
    df = df[cols]

    df = df.dropna(subset=['smiles']).reset_index(drop=True)
    try:
        refined_df = RefinementStage(df).complete_refinement(
            output_dir=OUT_DIR,
            get_report=False,
            n_cpu=1
        )

        if "smiles" in refined_df.columns:
            refined_df.rename(columns={"smiles": "SMILES"}, inplace=True)

        refined_df.to_csv(FINAL, index=False, encoding='utf-8')
        print(f"{FINAL}: {len(refined_df)} rows")

    except Exception as e:
        import traceback
        traceback.print_exc()   

if __name__ == "__main__":
    main()
