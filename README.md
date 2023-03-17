# Fully Functionalized Fragment Embeddings

Create a fragment embedding descriptor focused on the CRF fragment space. The descriptor is a lightweight 512-dimensional embedding capturing graph connectivity properties as well as physicochemical properties of the fragments.

## Installation

```bash
git clone https://github.com/ligand-discovery/fragment-embedding.git
cd fragment-embedding
python -m pip install -e .
```

## Usage

```python
from fragmentembedding import FragmentEmbedder

smiles_list = ["C#CCCC1(CCC(=O)N2CCN(Cc3cccc(F)c3)C(=O)C2)N=N1", "C#CCCC1(CCC(=O)N2CCN(Cc3cccc(C(F)(F)F)c3)CC2)N=N1", "C#CCCC1(CCC(=O)N2CCN3CC(F)(F)C[C@H]3C2)N=N1"]

fe = FragmentEmbedder()
X = fe.transform(smiles_list)

print(X.shape)
```

## About

This project was performed at [Georg Winter Lab](https://www.winter-lab.com/), based at [CeMM](https://cemm.at), Vienna.
