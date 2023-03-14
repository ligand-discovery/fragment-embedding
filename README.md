# Fully Functionalized Fragment Embeddings

Create a fragment embedding descriptor focused on the CRF fragment space. The descriptor is a lightweight 512-dimensional embedding capturing graph connectivity properties as well as physicochemical properties of the fragments.

## Installation

```bash
git clone https://github.com/ligand-discovery/fragment-descriptor.git
cd fragment-descriptor
python -m pip install -e .
```

## Usage

```python
from fragmentembedding import FragmentEmbedding

smiles_list = ["C#CCCC1(CCC(=O)N2CCN(Cc3cccc(F)c3)C(=O)C2)N=N1", "C#CCCC1(CCC(=O)N2CCN(Cc3cccc(C(F)(F)F)c3)CC2)N=N1", "C#CCCC1(CCC(=O)N2CCN3CC(F)(F)C[C@H]3C2)N=N1"]

fe = FragmentEmbedding()
X = fe.transform(smiles_list)
```
