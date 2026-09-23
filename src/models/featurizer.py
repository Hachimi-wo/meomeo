import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from config import MORGAN_RADIUS, MORGAN_BITS

_DESC_NAMES = None
_N_DESCRIPTORS = None

def _init_descriptors():
    global _DESC_NAMES, _N_DESCRIPTORS
    if _N_DESCRIPTORS is not None:
        return
    
    test_mol = Chem.MolFromSmiles("C")
    desc_dict = Descriptors.CalcMolDescriptors(test_mol)
    _DESC_NAMES = sorted(desc_dict.keys())
    _N_DESCRIPTORS = len(_DESC_NAMES)
    print(f"  [Featurizer] RDKit: {_N_DESCRIPTORS} 2D descriptors")


def get_n_descriptors():
    _init_descriptors()
    return _N_DESCRIPTORS


def compute_descriptors(mol):
    _init_descriptors()
    
    desc_dict = Descriptors.CalcMolDescriptors(mol)
    vals = []
    for name in _DESC_NAMES:
        v = desc_dict.get(name, 0.0)
        if v is None or np.isnan(v) or np.isinf(v):
            v = 0.0
        vals.append(float(v))
    return np.array(vals, dtype=np.float64)


def _count_fp_to_array(count_fp, n_bits):
    vec = np.zeros(n_bits, dtype=np.float32)
    for bit_id, count in count_fp.GetNonzeroElements().items():
        vec[bit_id % n_bits] = float(count)
    return vec


def featurize(smiles_list):
    _init_descriptors()
    
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=MORGAN_RADIUS, fpSize=MORGAN_BITS
    )

    feats, valid_idx = [], []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        parts = []

        parts.append(_count_fp_to_array(morgan_gen.GetCountFingerprint(mol), MORGAN_BITS))

        parts.append(compute_descriptors(mol))

        vec = np.concatenate(parts)
        feats.append(vec)
        valid_idx.append(i)

    X = np.array(feats, dtype=np.float32)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    dim_parts = [f"Morgan={MORGAN_BITS}"]
    expected = MORGAN_BITS
    dim_parts.append(f"Desc={_N_DESCRIPTORS}")
    expected += _N_DESCRIPTORS

    assert X.shape[1] == expected, \
        f"Dimension mismatch: {X.shape[1]} != {expected}"
    
    print(f"  [Featurizer] Feature vector: {X.shape[1]} dims "
          f"({' + '.join(dim_parts)})")
    return X, valid_idx
