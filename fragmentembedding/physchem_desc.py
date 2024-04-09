import numpy as np
import joblib
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.feature_selection import VarianceThreshold
from rdkit import Chem
import pandas as pd
from rdkit.Chem import Descriptors
from tqdm import tqdm


MAX_NA = 0.2


class NanFilter(object):
    def __init__(self):
        self._name = "nan_filter"

    def fit(self, X):
        max_na = int((1 - MAX_NA) * X.shape[0])
        idxs = []
        for j in range(X.shape[1]):
            c = np.sum(np.isnan(X[:, j]))
            if c > max_na:
                continue
            else:
                idxs += [j]
        self.col_idxs = idxs

    def transform(self, X):
        return X[:, self.col_idxs]

    def save(self, file_name):
        joblib.dump(self, file_name)

    def load(self, file_name):
        return joblib.load(file_name)


class Imputer(object):
    def __init__(self):
        self._name = "imputer"
        self._fallback = 0

    def fit(self, X):
        ms = []
        for j in range(X.shape[1]):
            vals = X[:, j]
            mask = ~np.isnan(vals)
            vals = vals[mask]
            if len(vals) == 0:
                m = self._fallback
            else:
                m = np.median(vals)
            ms += [m]
        self.impute_values = np.array(ms)

    def transform(self, X):
        for j in range(X.shape[1]):
            mask = np.isnan(X[:, j])
            X[mask, j] = self.impute_values[j]
        return X

    def save(self, file_name):
        joblib.dump(self, file_name)

    def load(self, file_name):
        return joblib.load(file_name)


class VarianceFilter(object):
    def __init__(self):
        self._name = "variance_filter"

    def fit(self, X):
        self.sel = VarianceThreshold()
        self.sel.fit(X)
        self.col_idxs = self.sel.transform([[i for i in range(X.shape[1])]]).ravel()

    def transform(self, X):
        return self.sel.transform(X)

    def save(self, file_name):
        joblib.dump(self, file_name)

    def load(self, file_name):
        return joblib.load(file_name)


def physchem_featurizer(smiles_list):
    R = []
    for smiles in tqdm(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        descriptors = []
        for _, descr_calc_fn in Descriptors._descList:
            descriptors.append(descr_calc_fn(mol))
        R += [np.array(descriptors)]
    return np.array(R)


def physchem_featurizer_as_dataframe(smiles_list):
    R = []
    for smiles in tqdm(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        descriptors = []
        for _, descr_calc_fn in Descriptors._descList:
            descriptors.append(descr_calc_fn(mol))
        R += [np.array(descriptors)]
    return pd.DataFrame(np.array(R), columns=[x[0] for x in Descriptors._descList])


class PhyschemDescriptor(object):
    def __init__(self, discretize=True):
        self.nan_filter = NanFilter()
        self.imputer = Imputer()
        self.variance_filter = VarianceFilter()
        self.discretizer = KBinsDiscretizer(
            n_bins=5, encode="ordinal", strategy="quantile"
        )
        self.discretize = discretize

    def fit(self, smiles):
        R = physchem_featurizer(smiles)
        X = np.array(R, dtype=np.float32)
        self.nan_filter.fit(X)
        X = self.nan_filter.transform(X)
        self.imputer.fit(X)
        X = self.imputer.transform(X)
        self.variance_filter.fit(X)
        X = self.variance_filter.transform(X)
        self.discretizer.fit(X)

    def transform(self, smiles):
        df = physchem_featurizer_as_dataframe(smiles)
        X = np.array(df, dtype=np.float32)
        X = self.nan_filter.transform(X)
        X = self.imputer.transform(X)
        X = self.variance_filter.transform(X)
        X = self.discretizer.transform(X)
        return np.array(X, dtype=int)


class PhyschemDescriptorWithFeatures(object):
    def __init__(self, discretize=True):
        self.nan_filter = NanFilter()
        self.imputer = Imputer()
        self.variance_filter = VarianceFilter()
        self.discretizer = KBinsDiscretizer(
            n_bins=5, encode="ordinal", strategy="quantile"
        )
        self.discretize = discretize

    def fit(self, smiles):
        df = physchem_featurizer_as_dataframe(smiles)
        X = np.array(df, dtype=np.float32)
        self.nan_filter.fit(X)
        X = self.nan_filter.transform(X)
        self.imputer.fit(X)
        X = self.imputer.transform(X)
        self.variance_filter.fit(X)
        X = self.variance_filter.transform(X)
        if self.discretize:
            self.discretizer.fit(X)
        col_idxs = self.variance_filter.col_idxs
        feature_names = list(df.columns)
        self.feature_names = [feature_names[i] for i in col_idxs]

    def transform(self, smiles):
        df = physchem_featurizer_as_dataframe(smiles)
        X = np.array(df, dtype=np.float32)
        X = self.nan_filter.transform(X)
        X = self.imputer.transform(X)
        X = self.variance_filter.transform(X)
        if self.discretize:
            X = self.discretizer.transform(X)
        return np.array(X, dtype=int)
