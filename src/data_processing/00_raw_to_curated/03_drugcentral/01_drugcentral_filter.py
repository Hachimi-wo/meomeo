# =============================================================================
# bỏ các INN đã xuất hiện trong ADReCS/Pneumotox.
# So khớp exact sau khi lower().strip().
# =============================================================================
import os
import pandas as pd

INTERMEDIATE_DIR = "../../../../data/01_intermediate/phase1_curation/drugcentral"
ADRECS = '../../../../data/02_curated_raw/adrecs_pneumotox_curated.csv'
DRUGCENTRAL = '../../../../data/00_raw_data/drugcentral/drugcentral.structures.smiles.tsv'
OUT = f'{INTERMEDIATE_DIR}/drugcentralonly.csv'

def main():
    pneu = pd.read_csv(ADRECS).dropna(subset=['Drug_Name'])
    dc = pd.read_csv(DRUGCENTRAL, sep='\t').dropna(subset=['INN'])

    # tên đã có bên ADReCS (chuẩn hoá lowercase + strip)
    pneu_set = set(pneu['Drug_Name'].astype(str).str.lower().str.strip())
    dc_names = dc['INN'].astype(str).str.lower().str.strip()

    # Giữ lại DrugCentral rows chưa xuất hiện trong ADReCS + pneumotox
    dc_only = dc[~dc_names.isin(pneu_set)].copy()

    os.makedirs(INTERMEDIATE_DIR, exist_ok=True)
    dc_only.to_csv(OUT, index=False, encoding='utf-8')
    print(f"{OUT}: {len(dc_only)} drugs")

if __name__ == "__main__":
    main()