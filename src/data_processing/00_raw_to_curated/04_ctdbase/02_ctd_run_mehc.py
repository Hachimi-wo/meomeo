# =============================================================================
#   Runs MEHC curation on CTD dataset to standardize SMILES.
# =============================================================================

import os, shutil
import pandas as pd
from mehc_curation.refinement import RefinementStage

INTERMEDIATE_DIR = "../../../../data/01_intermediate/phase1_curation/ctdbase"
INPUT = f"{INTERMEDIATE_DIR}/ctdbase_smiles.csv"
OUT_DIR = f"{INTERMEDIATE_DIR}/mehc"
FINAL = "../../../../data/02_curated_raw/ctd_mechanism_markers.csv"

def main():
    try:
        df = pd.read_csv(INPUT)        
        refined_df = RefinementStage(df).complete_refinement(
            output_dir=OUT_DIR, get_report=True
        )
        
    except Exception as e:
        print(f"[Error] An issue occurred during execution: {e}")

if __name__ == "__main__":
    main()