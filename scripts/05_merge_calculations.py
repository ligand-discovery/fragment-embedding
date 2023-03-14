import h5py
import os
from tqdm import tqdm
import numpy as np
from sklearn.model_selection import train_test_split

root = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.abspath(os.path.join(root, "..", "data"))

tsv_files = []
for l in os.listdir(data_folder):
    if "_split_" in l:
        tsv_files += [l]

with h5py.File(os.path.join(data_folder, "data_for_training.h5"), "w") as g:
    X = None
    for i in tqdm(range(len(tsv_files))):
        n = i + 1
        with h5py.File(
            os.path.join(data_folder, "grover", "fff_library_split_{0}.h5".format(n)),
            "r",
        ) as f:
            X_0 = f["V"][:]
        with h5py.File(
            os.path.join(data_folder, "mordred", "fff_library_split_{0}.h5".format(n)),
            "r",
        ) as f:
            X_1 = f["V"][:]
        X_ = np.hstack([X_0, X_1])
        if X is None:
            X = X_
        else:
            X = np.vstack([X, X_])
    idxs = [i for i in range(X.shape[0])]
    train_idxs, test_idxs = train_test_split(idxs, test_size=0.2)
    g.create_dataset("X", data=X.astype(np.float16))
    g.create_dataset("X_train", data=X[train_idxs].astype(np.float16))
    g.create_dataset("X_test", data=X[test_idxs].astype(np.float16))
    X = None
    for i in tqdm(range(len(tsv_files))):
        n = i + 1
        with h5py.File(
            os.path.join(data_folder, "morgan", "fff_library_split_{0}.h5".format(n)),
            "r",
        ) as f:
            X_0 = f["V"][:]
        with h5py.File(
            os.path.join(data_folder, "physchem", "fff_library_split_{0}.h5".format(n)),
            "r",
        ) as f:
            X_1 = f["V"][:]
        X_ = np.hstack([X_0, X_1])
        if X is None:
            X = X_
        else:
            X = np.vstack([X, X_])
    g.create_dataset("Y", data=X.astype(np.float16))
    g.create_dataset("Y_train", data=X[train_idxs].astype(np.float16))
    g.create_dataset("Y_test", data=X[test_idxs].astype(np.float16))
