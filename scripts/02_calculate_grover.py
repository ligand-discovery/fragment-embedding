import os
import csv
import pandas as pd
from tqdm import tqdm
import shutil
import h5py
import numpy as np

root = os.path.dirname(os.path.abspath(__file__))

import sys

sys.path.append(os.path.join(root, "..", "fragmentembedding"))

from grover_desc import GroverDescriptor

output_folder = os.path.join(root, "..", "data", "grover")
if os.path.exists(output_folder):
    shutil.rmtree(output_folder)
os.mkdir(output_folder)

desc = GroverDescriptor()
smiles = (
    pd.read_csv(os.path.join(root, "..", "data", "enamine_stock.csv"))[
        "smiles"
    ].tolist()
    + pd.read_csv(
        os.path.join(root, "..", "data", "fid2can_fff_all.tsv"), delimiter="\t"
    )["smiles"].tolist()
)

desc.fit(smiles)

for l in tqdm(os.listdir(os.path.join(root, "..", "data"))):
    if "_split_" in l:
        with open(os.path.join(root, "..", "data", l), "r") as f:
            reader = csv.reader(f)
            smiles = []
            for r in reader:
                smiles += r
        X = desc.transform(smiles)
        h5_file = os.path.join(output_folder, l.replace(".tsv", ".h5"))
        with h5py.File(h5_file, "w") as f:
            f.create_dataset("V", data=np.array(X))
