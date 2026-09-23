# =============================================================================
# Dedup theo Drug_Name, giữ đúng vị trí dòng gốc.
# nếu dòng CUỐI có "ADReCS" trong Pattern -> giữ dòng ĐẦU, ngược lại giữ CUỐI.
# Dùng csv module (không pandas) để bảo toàn thứ tự ghi ra.
# =============================================================================
import csv
from collections import defaultdict

IN  = "../../../data/01_intermediate/phase2_assembly/02_labeled_by_pattern.csv"
OUT = "../../../data/01_intermediate/phase2_assembly/03_deduplicated_by_drugname.csv"
KEY = "Drug_Name"

def main():
    groups = defaultdict(list)
    with open(IN, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
        for i, row in enumerate(reader):
            groups[row[KEY]].append((i, row))

    picked = []
    for rows in groups.values():
        first = min(rows, key=lambda x: x[0])
        last  = max(rows, key=lambda x: x[0])
        picked.append(first if "ADReCS" in last[1].get("Pattern", "") else last)

    picked.sort(key=lambda x: x[0])              # khôi phục thứ tự xuất hiện ban đầu

    with open(OUT, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=headers)
        w.writeheader()
        w.writerows(r for _, r in picked)

    total = sum(len(v) for v in groups.values())
    print(f"{OUT}: {total} -> {len(picked)} rows")

if __name__ == "__main__":
    main()