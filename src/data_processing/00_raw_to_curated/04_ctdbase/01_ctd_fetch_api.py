# =============================================================================
# Fetch CanonicalSMILES cho CTD drugs qua PubChem theo tên (thử nhiều variant).
# =============================================================================
import re, sys, time
from urllib.parse import quote
import requests, pandas as pd

INPUT_CSV  = "../../../../data/00_raw_data/ctdbase/need2find.csv"
OUTPUT_CSV = "../../../../data/01_intermediate/phase1_curation/ctdbase/ctdbase_smiles.csv"

SLEEP, TIMEOUT, RETRIES = 0.5, 15, 2

def clean_name(name):
    name = re.sub(r'\([^)]*\)', '', name)
    return re.sub(r'[,\s]+', ' ', name).strip()

def _try_name(n):
    if not n:
        return None
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           f"{quote(n)}/property/CanonicalSMILES/TXT")
    for att in range(RETRIES + 1):
        try:
            r = requests.get(url, timeout=TIMEOUT)
            if r.status_code == 200:
                s = r.text.strip()
                if s and s != "NotFound":
                    return s
            elif r.status_code == 404:
                return None
        except Exception:
            pass
        if att < RETRIES:
            time.sleep(1)
    return None

def get_smiles(name):
    variants = [name]
    if ',' in name:
        variants += [p.strip() for p in name.split(',')]
    cleaned = clean_name(name)
    if cleaned != name:
        variants.append(cleaned)
    variants = list(dict.fromkeys(v for v in variants if v))

    for v in variants:
        s = _try_name(v)
        if s:
            return s
    return None

def main():
    try:
        df = pd.read_csv(INPUT_CSV, encoding='utf-8')
    except FileNotFoundError:
        sys.exit(f"File not found: {INPUT_CSV}")

    names = df['Drug_Name'].dropna().unique().tolist()

    rows = []
    for name in names:
        rows.append({'Drug_Name': name, 'SMILES': get_smiles(name)})
        time.sleep(SLEEP)                       

    pd.DataFrame(rows).to_csv(OUTPUT_CSV, index=False, encoding='utf-8-sig')
    print(f"{OUTPUT_CSV}: {len(rows)} rows")

if __name__ == "__main__":
    main()