import pandas as pd
import csv
import os

df = pd.read_csv("../data/fid2can_fff_all.tsv", delimiter="\t")
smiles_0 = list(set(df["smiles"]))
df = pd.read_csv("../data/enamine_stock.csv", delimiter=",")
smiles_1 = list(set(df["smiles"]))
df = pd.read_csv("../data/enamine_real.csv", delimiter=",")
smiles_2 = list(set(df["smiles"]))

smiles = smiles_0 + smiles_1 + smiles_2

with open("../data/fff_library.tsv", "w") as f:
    writer = csv.writer(f)
    for smi in smiles:
        writer.writerow([smi])


def split(filehandler, delimiter=',', row_limit=10000,
          output_name_template='../data/fff_library_split_%s.tsv', output_path='.', keep_headers=False):
    reader = csv.reader(filehandler, delimiter=delimiter)
    current_piece = 1
    current_out_path = os.path.join(
        output_path,
        output_name_template % current_piece
    )
    current_out_writer = csv.writer(open(current_out_path, 'w'), delimiter=delimiter)
    current_limit = row_limit
    if keep_headers:
        headers = next(reader.next)
        current_out_writer.writerow(headers)
    for i, row in enumerate(reader):
        if i + 1 > current_limit:
            current_piece += 1
            current_limit = row_limit * current_piece
            current_out_path = os.path.join(
                output_path,
                output_name_template % current_piece
            )
            current_out_writer = csv.writer(open(current_out_path, 'w'), delimiter=delimiter)
            if keep_headers:
                current_out_writer.writerow(headers)
        current_out_writer.writerow(row)

split(open("../data/fff_library.tsv", "r"))