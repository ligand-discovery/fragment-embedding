import os
import sys
import pandas as pd
import numpy as np

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from fragment_embedder import FragmentEmbedder

fe = FragmentEmbedder()

df = pd.read_csv("../data/fid2can_fff_all.tsv", delimiter="\t")
smiles = df["smiles"].tolist()

X = fe.transform(smiles)
X = np.array(X)
print(X)
#np.save("../data/fid_fff_all_embedding", X)
