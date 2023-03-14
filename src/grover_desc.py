import os
import numpy as np
from ersilia import ErsiliaModel
from sklearn.preprocessing import RobustScaler


root = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.abspath(os.path.join(root, "..", "tools", "grover"))
MODEL = os.path.abspath(os.path.join(root, "..", "data", "checkpoints", "grover_large.pt"))


def grover_featurizer(smiles_list):
    with ErsiliaModel("grover-embedding") as mdl:
        X = mdl.predict(input=smiles_list, output="numpy")
    return X


class GroverDescriptor():

    def __init__(self):
        pass

    def fit(self, smiles):
        X = grover_featurizer(smiles)
        self.scaler = RobustScaler()
        self.scaler.fit(X)

    def transform(self, smiles):
        X = grover_featurizer(smiles)
        X = self.scaler.transform(X)
        return X
