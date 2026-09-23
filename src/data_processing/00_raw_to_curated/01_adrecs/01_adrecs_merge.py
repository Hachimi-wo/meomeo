# =============================================================================
#  Merges ADReCS and Pneumotox datasets, 
#keeping Pneumotox base and adding missing from ADReCS.
#Convert all drug names to lowercase and strip away salt designations (e.g., hydrochloride, sulfate) to merge the two databases.
#Assign the label "TO_BE_FOUND" to any drugs from Pneumotox that lack a code. 
# =============================================================================
import os, re, pandas as pd

INTERMEDIATE_DIR   = "../../../../data/01_intermediate/phase1_curation/adrecs"
BASE_EXT           = "../../../../data/00_raw_data"
ADRECS_FILE        = os.path.join(BASE_EXT, "adrecs", "ADReCS_Respiratory_Drugs.csv")
PNEUMOTOX_FILE     = os.path.join(BASE_EXT, "adrecs", "pneumotox_full_database.csv")
OUTPUT_FILE        = f"{INTERMEDIATE_DIR}/merged_database_v2.csv"
OUTPUT_ADRECS_ONLY = f"{INTERMEDIATE_DIR}/adrecs_only_drugs.csv"

REMOVE_TOKENS = [
    r'\bmonohydrate\b', r'\bdihydrate\b', r'\btrihydrate\b',
    r'\bhydrochloride\b', r'\bhydrobromide\b', r'\bmesylate\b',
    r'\bmaleate\b', r'\bfumarate\b', r'\bacetate\b', r'\bpropionate\b',
    r'\bsulfate\b', r'\bphosphate\b', r'\bcilexetil\b', r'\bexetil\b',
    r'\bproxetil\b', r'\bsodium\b', r'\bpotassium\b', r'\bcalcium\b',
    r'\bmagnesium\b', r'\bzinc\b', r'\b(?:i\.v\.|iv)\b',
]

def normalize_drug_name(name):
    if not isinstance(name, str):
        return ""
    name = name.lower().strip()
    for pat in REMOVE_TOKENS:
        name = re.sub(pat, '', name)
    return re.sub(r'\s+', ' ', name).strip()

def read_csv_any(path):
    try:
        return pd.read_csv(path, encoding='utf-8')
    except UnicodeDecodeError:
        return pd.read_csv(path, encoding='latin1')

def main():
    df_adrecs = read_csv_any(ADRECS_FILE)
    df_pneumo = read_csv_any(PNEUMOTOX_FILE)

    #  ADReCS 
    df_adrecs['norm'] = df_adrecs['DRUG_NAME'].apply(normalize_drug_name)
    df_adrecs = df_adrecs[df_adrecs['norm'] != '']
    grouped = df_adrecs.groupby('norm', as_index=False).agg({
        'DRUG_NAME':  lambda x: x.iloc[0],
        'PubChem_ID': lambda x: x.iloc[0] if len(x) > 0 else None,
        'ADR_TERM':   lambda x: '; '.join(sorted(set(x))),
    }).rename(columns={'DRUG_NAME': 'drug_name_rep', 'PubChem_ID': 'pubchem_id_rep'})
    adrecs_norms = set(grouped['norm'])

    #  Pneumotox 
    df_pneumo['norm'] = df_pneumo['Drug_Name'].apply(normalize_drug_name)
    df_pneumo = (df_pneumo[df_pneumo['norm'] != '']
                 .drop_duplicates(subset=['Drug_Name', 'Pattern']))
    pneumo_norms = set(df_pneumo['norm'])

    #  Group split 
    common_norms      = adrecs_norms & pneumo_norms
    adrecs_only_norms = adrecs_norms - pneumo_norms
    pneumo_only_norms = pneumo_norms - adrecs_norms

    #  ADReCS-only table 
    adrecs_only_df = grouped[grouped['norm'].isin(adrecs_only_norms)].copy()
    adrecs_only_df['Pattern'] = "ADReCS: " + adrecs_only_df['ADR_TERM']
    adrecs_only_df = (adrecs_only_df[['drug_name_rep', 'Pattern', 'pubchem_id_rep', 'norm']]
                      .rename(columns={'drug_name_rep': 'Drug_Name',
                                       'pubchem_id_rep': 'PubChem_ID'}))
    adrecs_only_df[['Drug_Name', 'norm', 'PubChem_ID']].to_csv(
        OUTPUT_ADRECS_ONLY, index=False, encoding='utf-8-sig')

    #  Merge
    merged_df = df_pneumo[['Drug_Name', 'Pattern', 'norm']].copy()
    pubchem_map = dict(zip(grouped['norm'], grouped['pubchem_id_rep']))
    merged_df['PubChem_ID'] = merged_df['norm'].map(pubchem_map)
    merged_df.loc[merged_df['norm'].isin(pneumo_only_norms), 'PubChem_ID'] = "TO_BE_FOUND"
    merged_df = pd.concat(
        [merged_df, adrecs_only_df[['Drug_Name', 'Pattern', 'norm', 'PubChem_ID']]],
        ignore_index=True)
    merged_df['Source'] = merged_df['norm'].map(
        lambda n: "Pneumotox_only" if n in pneumo_only_norms
                  else "ADReCS_only" if n in adrecs_only_norms else "Both")
    merged_df = merged_df[['Drug_Name', 'Pattern', 'PubChem_ID', 'Source', 'norm']]
    merged_df.to_csv(OUTPUT_FILE, index=False, encoding='utf-8-sig')

if __name__ == "__main__":
    main()
