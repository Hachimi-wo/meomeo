# =============================================================================
# 07_assemble_final_dataset.py
# Dataset cuối: positive (label có >=1 nhãn 1) + negative (score=0).
# Dedup theo SMILES keep first.
# =============================================================================
import pandas as pd

POOL = '../../../data/01_intermediate/phase2_assembly/05_merged_ctd_positive_pool.csv'
SCORED = '../../../data/01_intermediate/phase2_assembly/07_pubchem_scored_negatives.csv'
OUT = '../../../data/03_final_dataset/balanced_pulmonary_adr_dataset.csv'

LABEL_COLS = ['Label_ILD', 'Label_Edema', 'Label_Pleural', 'Label_Embolism']

def main():
    pool = pd.read_csv(POOL)
    positive_df = pool[pool[LABEL_COLS].any(axis=1)].copy()

    scored = pd.read_csv(SCORED)
    neg = scored[scored['score'] == 0]

    neg_rows = [{
        'Drug_Name':    r['Drug_Name'],
        'Pattern':      '',
        'PubChem_ID':   r['PubChem_ID'],
        'DrugBank_ID':  '',
        'KEGG_ID':      '',
        'SMILES':       r['SMILES'],
        'Label_ILD':      0,
        'Label_Edema':    0,
        'Label_Pleural':  0,
        'Label_Embolism': 0,
    } for _, r in neg.iterrows()]
    negative_df = pd.DataFrame(neg_rows)

    df = pd.concat([positive_df, negative_df], ignore_index=True)
    df = df[['Drug_Name', 'Pattern', 'PubChem_ID', 'DrugBank_ID', 'KEGG_ID',
             'SMILES'] + LABEL_COLS]
    df.drop_duplicates(subset=['SMILES'], keep='first', inplace=True)
    df.to_csv(OUT, index=False)
    print(f"positive={len(positive_df)} | negative={len(negative_df)} | total={len(df)} -> {OUT}")

if __name__ == "__main__":
    main()