# =============================================================================
# Đọc SIDER, lọc side-effect theo 4 pattern hô hấp, fetch SMILES  qua PubChem.
# =============================================================================
import csv, gzip, os, time
import requests

INTERMEDIATE_DIR = "../../../../data/01_intermediate/phase1_curation/sider"
EXTERNAL_DIR = "../../../../data/00_raw_data/sider"

TARGET_PATTERNS = {
    "Pleural effusion", "Pulmonary oedema",
    "Pulmonary embolism", "Interstitial pneumonia",
}
PATTERN_MAP = {p.lower(): p for p in TARGET_PATTERNS}   

def open_any(path, mode='rt'):
    if path.endswith('.gz'):
        return gzip.open(path, mode, encoding='utf-8')
    return open(path, mode, encoding='utf-8')

def match_pattern(side_effect):
    low = side_effect.lower()
    return {orig for pl, orig in PATTERN_MAP.items() if pl == low or pl in low}

def read_drug_names():
    path = f"{EXTERNAL_DIR}/drug_names.tsv"
    if not os.path.exists(path):
        return {}
    with open_any(path) as f:
        return {r[0].strip(): r[1].strip()
                for r in csv.reader(f, delimiter='\t') if len(r) >= 2}

def collect_patterns(path):
    out = {}
    with open_any(path) as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) < 6:
                continue
            matched = match_pattern(row[5].strip())
            if matched:
                out.setdefault(row[0].strip(), set()).update(matched)
    return out

def fetch_smiles(cid):
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
           f"cid/{cid}/property/CanonicalSMILES/TXT")
    try:
        r = requests.get(url, timeout=15)
        return r.text.strip() if r.status_code == 200 else "NA"
    except requests.exceptions.RequestException:
        return "NA"

def main():
    out = f"{INTERMEDIATE_DIR}/sider_filtered_nomehc.csv"
    
    # Load cache if exists
    cached_results = {}
    if os.path.exists(out):
        with open(out, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                cached_results[row["STITCH_ID"]] = row

    drug_names = read_drug_names()
    drug_patterns = collect_patterns(f"{EXTERNAL_DIR}/meddra_all_se.tsv")

    results = []
    for stitch_id, patterns in drug_patterns.items():
        if stitch_id in cached_results:
            results.append(cached_results[stitch_id])
            continue
            
        cid = stitch_id[3:] if stitch_id.startswith('CID') else stitch_id
        results.append({
            "STITCH_ID": stitch_id,
            "Drug_Name": drug_names.get(stitch_id, f"STITCH_{stitch_id}"),
            "SMILES": fetch_smiles(cid),
            "Respiratory_Patterns": "; ".join(sorted(patterns)),
        })
        time.sleep(1.5)     

    with open(out, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=["STITCH_ID", "Drug_Name",
                                          "SMILES", "Respiratory_Patterns"])
        w.writeheader()
        w.writerows(results)
    print(f"{out}: {len(results)} drugs")

if __name__ == "__main__":
    main()