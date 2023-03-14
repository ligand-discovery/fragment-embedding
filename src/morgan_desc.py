from rdkit import Chem, DataStructs
import numpy as np
from rdkit.Chem import rdFingerprintGenerator
from sklearn.feature_selection import VarianceThreshold


def get_ecfp_fingerprint(smiles_list):
    R = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        fingerprints_vect = rdFingerprintGenerator.GetCountFPs(
            [mol], fpType=rdFingerprintGenerator.MorganFP
        )[0]
        fingerprint = np.zeros((0,), np.float32)  # Generate target pointer to fill
        DataStructs.ConvertToNumpyArray(fingerprints_vect, fingerprint)
        R += [fingerprint]
    X = np.array(R, dtype=int)
    return X


class MorganFingerprint(object):

    def __init__(self):
        self.variance_filter = VarianceThreshold(threshold=0)

    def fit(self, smiles):
        X = get_ecfp_fingerprint(smiles)
        self.variance_filter.fit(X)

    def transform(self, smiles):
        X = get_ecfp_fingerprint(smiles)
        X = self.variance_filter.transform(X)
        return X
