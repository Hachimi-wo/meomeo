# =============================================================================
# - Resolve PubChem CID cho drug "TO_BE_FOUND" chưa có CID .
# - Fetch CanonicalSMILES cho tất cả CID hợp lệ.
# =============================================================================
import os, json, time, urllib.parse
import requests, pandas as pd
from tqdm import tqdm

INTERMEDIATE_DIR = "../../../../data/01_intermediate/phase1_curation/adrecs"
INPUT_FILE = f'{INTERMEDIATE_DIR}/merged_database_v2.csv'
OUTPUT_FILE = f'{INTERMEDIATE_DIR}/adrecs_pubchem_smiles.csv'
SMILES_CP = f'{INTERMEDIATE_DIR}/smiles_checkpoint.csv'
CID_CP = f'{INTERMEDIATE_DIR}/cid_checkpoint.csv'

BATCH, TIMEOUT, RETRIES, DELAY = 5, 90, 6, 1.0
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
HEADERS = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) '
                         'AppleWebKit/537.36 (KHTML, like Gecko) '
                         'Chrome/91.0.4472.124 Safari/537.36'}

def _json(text):
    try: return json.loads(text)
    except json.JSONDecodeError: return None

# Tìm CID từ tên thuốc
def cid_from_name(name):
    try:
        url = f"{PUBCHEM}/name/{urllib.parse.quote(name)}/cids/JSON"
        r = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
        if r.status_code == 200:
            cids = r.json().get('IdentifierList', {}).get('CID')
            if cids: return cids[0]
    except Exception:
        pass
    return None

# Fallback : fetch từng CID một
def smiles_one(cid):
    try:
        url = f"{PUBCHEM}/cid/{cid}/property/CanonicalSMILES/JSON"
        r = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
        if r.status_code == 200:
            props = (_json(r.text) or {}).get('PropertyTable', {}).get('Properties', [])
            if props: return props[0].get('CanonicalSMILES', '')
    except Exception:
        pass
    return ''

def smiles_batch(cids, done=None):
    done = done or set()
    mapping = {}
    for i in tqdm(range(0, len(cids), BATCH), desc="SMILES"):
        batch = cids[i:i+BATCH]
        if all(str(c) in done for c in batch):
            continue
        cid_str = ','.join(str(c) for c in batch)
        url = f"{PUBCHEM}/cid/{cid_str}/property/CanonicalSMILES,ConnectivitySMILES/JSON"
        ok = False
        for att in range(RETRIES):
            try:
                r = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
                if r.status_code == 200:
                    props = (_json(r.text) or {}).get('PropertyTable', {}).get('Properties')
                    if props:
                        for p in props:
                            mapping[str(p['CID'])] = p.get('CanonicalSMILES',
                                                            p.get('ConnectivitySMILES', ''))
                        ok = True
                        break
                    time.sleep(5 * 2**att)      
                elif r.status_code == 503:
                    time.sleep(10 * 2**att)
                elif r.status_code == 404:
                    break
                else:
                    time.sleep(5)
            except Exception:
                time.sleep(5)
        if not ok:                               
            for cid in batch:
                s = smiles_one(cid)
                if s: mapping[str(cid)] = s
                time.sleep(0.5)
        time.sleep(DELAY)
        pd.DataFrame(mapping.items(), columns=['CID', 'SMILES']).to_csv(
            SMILES_CP, index=False)
    return mapping

def main():
    df = pd.read_csv(INPUT_FILE, encoding='utf-8-sig')
    print(f"{len(df)} rows from {INPUT_FILE}")

    # resolve CID cho các TO_BE_FOUND
    cid_cp = {}
    if os.path.exists(CID_CP):
        cp = pd.read_csv(CID_CP)
        cid_cp = dict(zip(cp['Drug_Name'], cp['PubChem_ID']))

    todo = df[df['PubChem_ID'] == 'TO_BE_FOUND']['Drug_Name'].unique()
    resolved = {}
    if len(todo):
        print(f"Resolving {len(todo)} CID via name...")
        for name in tqdm(todo):
            if name in cid_cp:
                resolved[name] = cid_cp[name]
                continue
            resolved[name] = cid_from_name(name) or "NOT_FOUND"
            time.sleep(DELAY)
            cid_cp.update(resolved)
            pd.DataFrame(cid_cp.items(), columns=['Drug_Name', 'PubChem_ID']).to_csv(
                CID_CP, index=False)

    df['PubChem_ID'] = df.apply(
        lambda r: resolved.get(r['Drug_Name'], "NOT_FOUND")
        if r['PubChem_ID'] == 'TO_BE_FOUND' else r['PubChem_ID'], axis=1)

    # fetch SMILES cho tất cả CID hợp lệ 
    cids = []
    for c in df['PubChem_ID'].unique():
        try: cids.append(int(float(c)))
        except Exception: pass
    print(f"{len(cids)} unique valid CIDs")

    done, cached = set(), {}
    if os.path.exists(SMILES_CP):
        cp = pd.read_csv(SMILES_CP)
        cp['SMILES'] = cp['SMILES'].fillna('')
        cached = dict(zip(cp['CID'].astype(str), cp['SMILES']))
        done = set(cached.keys())
        print(f"Checkpoint: {len(done)} CIDs")

    mapping = smiles_batch(cids, done)
    mapping.update(cached)

    def lookup(x):
        if pd.isna(x): return ''
        s = str(x).replace('.0', '')
        return mapping.get(s, '') if s.isdigit() else ''
    df['SMILES'] = df['PubChem_ID'].apply(lookup)

    df.to_csv(OUTPUT_FILE, index=False, encoding='utf-8-sig')
    print(f"SMILES covered: {(df['SMILES'] != '').sum()}/{len(df)}")

if __name__ == "__main__":
    main()