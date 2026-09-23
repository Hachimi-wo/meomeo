# =============================================================================
# Thêm DrugCentral (INN + SMILES) làm negative candidates (label=0).
# =============================================================================
import pandas as pd

IN      = '../../../data/01_intermediate/phase2_assembly/03_deduplicated_by_drugname.csv'
DC      = '../../../data/02_curated_raw/drugcentral_smiles.csv'
OUT     = '../../../data/01_intermediate/phase2_assembly/04_augmented_drugcentral.csv'

def main():
    df_orig = pd.read_csv(IN, encoding='utf-8')
    df_new  = pd.read_csv(DC, encoding='utf-8')
    if 'Unnamed: 0' in df_new.columns:
        df_new = df_new.drop(columns=['Unnamed: 0'])

    positive_drugs = set(df_orig['Drug_Name'].str.lower())
    df_new = df_new[~df_new['INN'].str.lower().isin(positive_drugs)]

    df_new_labels = pd.DataFrame({
        'Drug_Name':   df_new['INN'],
        'Pattern':     '',
        'PubChem_ID':  None,
        'DrugBank_ID': None,
        'KEGG_ID':     None,
        'SMILES':      df_new['SMILES'],
        'Label_ILD':      0,
        'Label_Edema':    0,
        'Label_Embolism': 0,
        'Label_Pleural':  0,
    })

    pd.concat([df_orig, df_new_labels], ignore_index=True) \
      .to_csv(OUT, index=False, encoding='utf-8')
    print(f"{OUT}: {len(df_orig)} + {len(df_new_labels)}")

if __name__ == "__main__":
    main()
