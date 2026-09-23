# =============================================================================
# Gộp ADReCS/PNEUMOTOX + SIDER. SIDER thiếu PubChem_ID/DrugBank_ID/KEGG_ID -> NaN.
# Output giữ đúng cột của ADReCS (df_final.columns).
# =============================================================================
import os
import pandas as pd

SIDER = "../../../data/02_curated_raw/sider_respiratory_adrs.csv"
ADRECS = "../../../data/02_curated_raw/adrecs_pneumotox_curated.csv"
OUT = "../../../data/01_intermediate/phase2_assembly/01_merged_adrecs_sider.csv"

def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    df_adrecs = pd.read_csv(ADRECS)
    df_sider = pd.read_csv(SIDER)[['Drug_Name', 'Respiratory_Patterns', 'SMILES']] \
                                 .rename(columns={'Respiratory_Patterns': 'Pattern'})

    df = pd.concat([df_adrecs, df_sider], ignore_index=True)
    df = df[df_adrecs.columns.tolist()]     
    df.to_csv(OUT, index=False)
    print(f"{OUT}: {len(df_adrecs)} + {len(df_sider)} -> {len(df)}")

if __name__ == "__main__":
    main()