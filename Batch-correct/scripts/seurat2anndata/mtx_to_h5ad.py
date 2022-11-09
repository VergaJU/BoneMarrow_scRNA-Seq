#!/home/jacopo/miniconda3/envs/scanpy_env/bin/python3

import scanpy as sc
import pandas as pd
import sys

adata = sc.read_10x_mtx("./temp/") # load data (10x matrix obtained with seurat_to_mtx.R)

meta = pd.read_csv("./temp/metadata.tsv", sep = "\t") # load metadata

# add metadata to adata object
for column in meta:
    adata.obs[column] = meta[column]


# set new filename with h5ad suffix
filename = str(sys.argv[1])[:-3]+"h5ad"

# safe h5ad file
adata.write(filename)

