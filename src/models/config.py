import os

CSV_PATH = "balanced_pulmonary_adr_dataset.csv"
SMILES_COL = "SMILES"
LABEL_COLS = ["Label_ILD", "Label_Edema", "Label_Pleural", "Label_Embolism"]

FRAC_TRAIN = 0.70
FRAC_VAL = 0.15
FRAC_TEST = 0.15

CV_FOLDS = 5

MORGAN_RADIUS, MORGAN_BITS = 2, 2048



VARIANCE_THRESHOLD = 0.01
SELECT_K_BEST = 300

N_ESTIMATORS = 200

OUT_DIR = "outputs"
os.makedirs(OUT_DIR, exist_ok=True)

def get_feature_groups():
    from featurizer import get_n_descriptors
    n_desc = get_n_descriptors()
    
    offset = 0
    groups = {}
    
    groups['morgan'] = (offset, offset + MORGAN_BITS)
    offset += MORGAN_BITS
    
    groups['desc'] = (offset, offset + n_desc)
    offset += n_desc
        
    return groups