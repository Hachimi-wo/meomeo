# =============================================================================
# Gộp CTD markers (Pattern = "ILD"/"pulmonary edema"/...) vào positive pool.
# =============================================================================
import pandas as pd

RUN2 = "../../../data/01_intermediate/phase2_assembly/04_augmented_drugcentral.csv"
CTD  = "../../../data/02_curated_raw/ctd_mechanism_markers.csv"
OUT  = "../../../data/01_intermediate/phase2_assembly/05_merged_ctd_positive_pool.csv"

LABEL_COLS = ["Label_ILD", "Label_Edema", "Label_Pleural", "Label_Embolism"]

def ctd_labels(pattern):
    p = str(pattern).strip().casefold()
    return {
        "ild":                {"Label_ILD": 1, "Label_Edema": 0, "Label_Pleural": 0, "Label_Embolism": 0},
        "pulmonary edema":    {"Label_ILD": 0, "Label_Edema": 1, "Label_Pleural": 0, "Label_Embolism": 0},
        "pleural effusion":   {"Label_ILD": 0, "Label_Edema": 0, "Label_Pleural": 1, "Label_Embolism": 0},
        "pulmonary embolism": {"Label_ILD": 0, "Label_Edema": 0, "Label_Pleural": 0, "Label_Embolism": 1},
    }.get(p, {c: 0 for c in LABEL_COLS})

def main():
    df_run2 = pd.read_csv(RUN2, encoding='utf-8')
    df_ctd  = pd.read_csv(CTD,  encoding='utf-8')

    label_df = pd.DataFrame(df_ctd["Pattern"].apply(ctd_labels).tolist())
    df_ctd_labeled = pd.concat([df_ctd.reset_index(drop=True), label_df], axis=1)

    for c in LABEL_COLS:
        if c not in df_run2.columns:
            df_run2[c] = 0

    df = pd.concat([df_run2, df_ctd_labeled], ignore_index=True, sort=False)
    df.to_csv(OUT, index=False, encoding='utf-8')
    print(f"{OUT}: {len(df_run2)} + {len(df_ctd)} -> {len(df)}")

if __name__ == "__main__":
    main()