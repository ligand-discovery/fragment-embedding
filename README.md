# Fully Functionalized Fragment Embeddings

This repository contains a lightweight, pre-trained neural network to obtain **512-dimensional embeddings** (vectors) specific to fully functionalized fragments.

The model has been trained on the Fully Functionalized Fragment subsets of the [Enamine Stock](https://enamine.net/compound-libraries) and [Enamine REAL](https://enamine.net/compound-collections/real-compounds/real-database-subsets) databases, adding up to almost 300,000 compounds.

These embeddings are aimed at capturing graph connectivity properties as well as physicochemical properties of the fully functionalized fragments screened in the [Ligand Discovery](https://ligand-discovery.ai) project.

## Installation

```bash
git clone https://github.com/ligand-discovery/fragment-embedding.git
cd fragment-embedding
python -m pip install -e .
```

Please note that installation should be done in development `-e` mode to ensure that model files (`.joblib`) are accessible to the package.

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

This work is done by the [Georg Winter Lab](https://www.winter-lab.com/) at [CeMM](https://cemm.at), Vienna.
