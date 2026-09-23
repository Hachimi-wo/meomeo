# =============================================================================
# Bổ sung SMILES còn thiếu qua ChEMBL (tra bằng DrugBank_ID).
# =============================================================================
import time
import pandas as pd
import requests
from tqdm import tqdm

INTERMEDIATE_DIR = "../../../../data/01_intermediate/phase1_curation/adrecs"
INPUT_FILE = f'{INTERMEDIATE_DIR}/adrecs_pubchem_smiles.csv'
OUTPUT_FILE = f'{INTERMEDIATE_DIR}/adrecs_chembl_smiles.csv'

TIMEOUT, RETRIES, DELAY = 30, 3, 0.3
HEADERS = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) '
                         'AppleWebKit/537.36 (KHTML, like Gecko) '
                         'Chrome/91.0.4472.124 Safari/537.36'}

def smiles_from_chembl(db_id):
    if not db_id or not db_id.startswith('DB'):
        return ''
    url = (f"https://www.ebi.ac.uk/chembl/api/data/molecule"
           f"?molecule_dictionary__drugbank_id={db_id}&format=json")
    for att in range(RETRIES):
        try:
            r = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
            if r.status_code == 200:
                mols = r.json().get('molecules', [])
                if mols:
                    return mols[0].get('molecule_structures', {}).get('canonical_smiles', '')
                return ''
            if r.status_code == 404:
                return ''
            time.sleep(1)
        except requests.exceptions.Timeout:
            time.sleep(2 ** att)
        except Exception:
            return ''
    return ''

def main():
    try:
        df = pd.read_csv(INPUT_FILE, encoding='utf-8-sig')
    except FileNotFoundError:
        return

    if 'SMILES' not in df.columns:
        return

    todo = df[df['SMILES'].isna() | (df['SMILES'] == '')].copy()

    if 'DrugBank_ID' not in df.columns or len(todo) == 0:
        df.to_csv(OUTPUT_FILE, index=False, encoding='utf-8-sig')
        return

    n = 0
    for idx, row in tqdm(todo.iterrows(), total=len(todo), desc="ChEMBL"):
        s = smiles_from_chembl(str(row['DrugBank_ID']).strip())
        if s:
            df.at[idx, 'SMILES'] = s
            n += 1
        time.sleep(DELAY)

    df.to_csv(OUTPUT_FILE, index=False, encoding='utf-8-sig')
    print(f"+{n} SMILES | total {(df['SMILES'] != '').sum()}/{len(df)}")

if __name__ == "__main__":
    main()