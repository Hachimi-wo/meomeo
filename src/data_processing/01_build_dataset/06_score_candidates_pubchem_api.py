# =============================================================================
# Chấm điểm candidate negatives bằng PubChem PUG-View + Assay Summary.
# Score = số lần keyword xuất hiện trong các section "adverse/toxicity/...".
# Có checkpoint (checkpoint.txt + results_partial.csv trong cwd) để resume.
# =============================================================================
import os, re, json, time
from urllib.parse import quote
import requests, pandas as pd
from tqdm import tqdm

IN  = "../../../data/01_intermediate/phase2_assembly/06_candidate_negatives_raw.csv"
OUT = "../../../data/01_intermediate/phase2_assembly/07_pubchem_scored_negatives.csv"
CP_FILE, CP_PARTIAL, CP_EVERY = "checkpoint.txt", "results_partial.csv", 50
SLEEP, TIMEOUT = 0.3, 15

KEYWORDS = [
    "interstitial lung disease", "ild", "pulmonary fibrosis", "pneumonitis",
    "interstitial pneumonia", "fibrosing alveolitis",
    "usual interstitial pneumonia", "desquamative interstitial pneumonia",
    "nonspecific interstitial pneumonia", "cryptogenic organizing pneumonia",
    "hypersensitivity pneumonitis", "sarcoidosis",
    "idiopathic pulmonary fibrosis", "ipf",
    "pleural effusion", "pleurisy", "pleuritis", "pneumothorax",
    "pleural thickening", "pleural plaque", "mesothelioma",
    "pulmonary edema", "pulmonary oedema", "acute pulmonary edema",
    "noncardiogenic pulmonary edema", "fluid accumulation", "hydrothorax",
    "lung edema", "capillary leak",
    "pulmonary embolism", "thromboembolism", "venous thromboembolism",
    "pulmonary thromboembolism", "deep vein thrombosis", "dvt", "pe",
]
ADVERSE_TERMS = ["adverse", "side effect", "toxicity", "safety", "hazard", "warning",
                 "drug induced", "pulmonary toxicity", "respiratory toxicity", "lung toxicity"]

def _clean(s):
    if pd.isna(s): return None
    s = str(s).strip()
    return None if s in ("", "NAN", "None") else s

def cid_from_pubchem_id(pid):
    s = _clean(pid)
    if s is None or s == "NOT_FOUND": return None
    try:    return int(float(s))
    except ValueError: return None

def _query_cid(name):
    r = requests.get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{quote(name)}/cids/JSON", timeout=TIMEOUT)
    if r.status_code == 200:
        cids = r.json().get("IdentifierList", {}).get("CID", [])
        if cids: return cids[0]
    return None

def cid_from_name(name):
    """Filter ký tự đặc biệt trước; nếu rỗng thì dùng nguyên tên."""
    n = _clean(name)
    if n is None: return None
    filtered = re.sub(r'[^a-zA-Z0-9\s\-]', '', n).strip() or n
    try:    return _query_cid(filtered)
    except Exception as e:
        print(f"      Lỗi API name: {e}"); return None

def cid_from_name_fallback(name):
    n = _clean(name)
    if n is None: return None
    try:    return _query_cid(n)
    except Exception: return None

def pug_view(cid):
    try:
        r = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON",
                         timeout=TIMEOUT)
        return r.json() if r.status_code == 200 else None
    except Exception: return None

def assay_summary(cid):
    try:
        r = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON",
                         timeout=TIMEOUT)
        return r.json() if r.status_code == 200 else None
    except Exception: return None

def _count(text, kws):
    return sum(text.count(k.lower()) for k in kws)

def score_pug(data):
    if not data: return 0
    def walk(sec):
        cnt = 0
        if any(t in sec.get("TOCHeading", "").lower() for t in ADVERSE_TERMS):
            for info in sec.get("Information", []):
                cnt += _count(json.dumps(info).lower(), KEYWORDS)
        for sub in sec.get("Section", []):
            cnt += walk(sub)
        return cnt
    secs = data.get("Record", {}).get("Section", []) or data.get("Section", [])
    return sum(walk(s) for s in secs)

def score_assay(data):
    if not data: return 0
    assays = data.get("Assays", []) or (data if isinstance(data, list) else [])
    return sum(_count(json.dumps(a).lower(), KEYWORDS)
               for a in assays
               if any(t in json.dumps(a).lower() for t in ADVERSE_TERMS))

def load_cp():
    if os.path.exists(CP_FILE):
        return int(open(CP_FILE).read().strip())
    return 0

def save_cp(i): open(CP_FILE, "w").write(str(i))

def load_partial():
    return pd.read_csv(CP_PARTIAL).to_dict('records') if os.path.exists(CP_PARTIAL) else []

def save_partial(rows): pd.DataFrame(rows).to_csv(CP_PARTIAL, index=False)

def main():
    df = pd.read_csv(IN)
    start = load_cp()
    results = load_partial()

    for idx in tqdm(range(start, len(df)), initial=start, total=len(df), desc="Xử lý"):
        row = df.iloc[idx]
        name, smiles, pid = row.get("Drug_Name", ""), row.get("SMILES", ""), row.get("PubChem_ID", "")

        # 1) thử lấy CID từ PubChem_ID -> tên (filtered) -> tên (fallback)
        cid, src = cid_from_pubchem_id(pid), "PubChem_ID"
        if not cid:
            cid, src = cid_from_name(name), "Drug_Name"
            if not cid:
                cid, src = cid_from_name_fallback(name), "Drug_Name (fallback)"
        if not cid:
            results.append({"Drug_Name": name, "SMILES": smiles, "PubChem_ID": pid,
                            "CID": None, "CID_Source": "none", "score": 0,
                            "status": "no_cid", "PUG_View": 0, "Assay_Summary": 0})
            time.sleep(SLEEP)
            if (idx + 1) % CP_EVERY == 0:
                save_cp(idx + 1); save_partial(results)
            continue

        # 2) chấm điểm qua PUG-View + Assay Summary
        pv, ascore = pug_view(cid), 0
        spv = score_pug(pv) if pv else 0
        sa  = assay_summary(cid)
        sas = score_assay(sa) if sa else 0
        total = spv + sas

        results.append({"Drug_Name": name, "SMILES": smiles, "PubChem_ID": pid,
                        "CID": cid, "CID_Source": src, "score": total,
                        "status": "ok" if (pv or sa) else "no_data",
                        "PUG_View": spv, "Assay_Summary": sas})

        if (idx + 1) % CP_EVERY == 0:
            save_cp(idx + 1); save_partial(results)
        time.sleep(SLEEP)

    save_cp(len(df)); save_partial(results)

    df_res = pd.DataFrame(results)
    df_res.to_csv(OUT, index=False)
    print(f"{OUT}: {len(df_res)} rows | có CID: {df_res['CID'].notna().sum()} | "
          f"score min/max/mean = {df_res['score'].min()}/"
          f"{df_res['score'].max()}/{df_res['score'].mean():.2f}")

    for f in (CP_FILE, CP_PARTIAL):
        if os.path.exists(f): os.remove(f)

if __name__ == "__main__":
    main()