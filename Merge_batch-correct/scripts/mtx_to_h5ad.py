#!/home/jacopo/miniconda3/envs/scanpy_env/bin/python3

import scanpy as sc
import pandas as pd
import sys

adata = sc.read_10x_mtx("./temp/") # load data (10x matrix obtained with seurat_to_mtx.R)

meta = pd.read_csv("./temp/metadata.tsv", sep = "\t") # load metadata

# add metadata to adata object
adata.obs["nCount_RNA"] = meta["nCount_RNA"]
adata.obs["nFeature_RNA"] = meta["nFeature_RNA"]
adata.obs["percent_mt"] = meta["percent_mt"]

# set new filename with h5ad suffix
filename = str(sys.argv[1])[:-3]+"h5ad"

# safe h5ad file
adata.write(filename)

