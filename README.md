# Fully Functionalized Fragment Embeddings

This repository contains a lightweight, pre-trained neural network to obtain a 512-dimensional embedding (vector) specific to fully functionalized fragments. The model has been trained on the Enamine REAL and Enamine Stock FFF libraries. It is aimed at capturing graph connectivity properties as well as physicochemical properties of the CRF fragments.

## Installation

```bash
git clone https://github.com/ligand-discovery/fragment-embedding.git
cd fragment-embedding
python -m pip install -e .
```

## Usage

Fragment embeddings can be produced with a simple Python API as follows:

```python
from fragmentembedding import FragmentEmbedder

smiles_list = [
        "C#CCCC1(CCC(=O)N2CCN(Cc3cccc(F)c3)C(=O)C2)N=N1",
        "C#CCCC1(CCC(=O)N2CCN(Cc3cccc(C(F)(F)F)c3)CC2)N=N1",
        "C#CCCC1(CCC(=O)N2CCN3CC(F)(F)C[C@H]3C2)N=N1"
    ]

fe = FragmentEmbedder()
X = fe.transform(smiles_list)

print(X.shape)
```

The expected shape of `X` is `(3, 512)`.

## About

This project was performed at [Georg Winter Lab](https://www.winter-lab.com/), based at [CeMM](https://cemm.at), Vienna.
