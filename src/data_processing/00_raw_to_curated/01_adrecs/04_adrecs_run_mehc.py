# =============================================================================
#   Runs MEHC curation pipeline on dataset to standardize SMILES.
# =============================================================================

import pandas as pd
INTERMEDIATE_DIR = "../../../../data/01_intermediate/phase1_curation/adrecs"

from mehc_curation.refinement import RefinementStage

def run_curation_pipeline():
    file_path = f"{INTERMEDIATE_DIR}/adrecs_chembl_smiles.csv"
    df = pd.read_csv(file_path)
    
    if "SMILES" in df.columns:
        df.rename(columns={"SMILES": "smiles"}, inplace=True)
        
    cols = ['smiles'] + [col for col in df.columns if col != 'smiles']
    df = df[cols]
    
    df = df.dropna(subset=['smiles']).reset_index(drop=True)
    
    refiner = RefinementStage(df)
    
    try:
        refined_df = refiner.complete_refinement(
            output_dir=f"{INTERMEDIATE_DIR}/mehc",
            get_report=True,
            n_cpu=1
        )        
    except Exception as e:
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    run_curation_pipeline()