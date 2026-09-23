# =============================================================================
# 02_annotate_labels_from_pattern.py
# Gán 4 nhãn từ cột Pattern theo 3 tầng: code (I.b, II.a...) -> ADReCS clause -> text.
# LƯU Ý: giữ nguyên so khớp 'ILD' in text (không bao giờ match vì text đã lower).
# =============================================================================
import re
import pandas as pd

IN  = '../../../data/01_intermediate/phase2_assembly/01_merged_adrecs_sider.csv'
OUT = '../../../data/01_intermediate/phase2_assembly/02_labeled_by_pattern.csv'

CODE_RE = re.compile(r'[A-Z]+[a-z]*\.[a-z]+')            
ADR_RE  = re.compile(r'ADReCS:\s*(.*?)(?=\n|$)', re.IGNORECASE)

def extract_codes(pattern):
    return [] if pd.isna(pattern) else CODE_RE.findall(pattern)

def extract_clauses(pattern):
    if pd.isna(pattern): return []
    return [m.strip() for m in ADR_RE.findall(pattern) if m.strip()]

def _hit(codes, clauses_low, pat_low, code_prefixes, text_keys):
    # 1) khớp theo code
    for c in codes:
        if c.startswith(tuple(code_prefixes)) or c in code_prefixes:
            return 1
    # 2) khớp trong ADReCS clause (đã lower)
    for cl in clauses_low:
        if any(k in cl for k in text_keys):
            return 1
    # 3) fallback khớp cả pattern gốc
    if any(k in pat_low for k in text_keys):
        return 1
    return 0

def assign_labels(pattern):
    pat_low = str(pattern).lower() if not pd.isna(pattern) else ''
    codes   = extract_codes(pattern)
    clauses = [c.lower() for c in extract_clauses(pattern)]

    ild      = _hit(codes, clauses, pat_low, ['I.b', 'I.g', 'g'], ['interstitial', 'fibrosis', 'ILD'])
    edema    = _hit(codes, clauses, pat_low, ['II.a', 'II.d'],
                    ['pulmonary oedema', 'pulmonary edema',
                     'noncardiogenic pulmonary edema', 'noncardiogenic pulmonary oedema'])
    pleural  = _hit(codes, clauses, pat_low, ['V.a'], ['pleural'])
    embolism = _hit(codes, clauses, pat_low, ['VI.a', 'VI.b', 'VI.d'], ['embolism'])
    return ild, edema, pleural, embolism

def main():
    df = pd.read_csv(IN, encoding='utf-8')

    labels = df['Pattern'].apply(assign_labels)
    df['Label_ILD']      = [x[0] for x in labels]
    df['Label_Edema']    = [x[1] for x in labels]
    df['Label_Pleural']  = [x[2] for x in labels]
    df['Label_Embolism'] = [x[3] for x in labels]

    df.to_csv(OUT, index=False, encoding='utf-8')
    print(f"{OUT}: {len(df)} rows")

if __name__ == "__main__":
    main()