import pandas as pd
import numpy as np
import scanpy as sc
import pyswarms as ps
import scib


# Read file
adata = sc.read("./SMM.h5ad")


# filter genes expressed in less than 3 cells
sc.pp.filter_genes(adata, min_cells = 3)


# normalize data
sc.pp.normalize_total(adata, target_sum=1e4)

adata_norm = adata.X

# convert to dense matrix
adata_norm = adata_norm.toarray()

# random matrix with same dimension:
rows = adata_norm.shape[0]
cols = adata_norm.shape[1]

random_matrix = np.random.rand(rows,cols)

# add matrices
new_matrix = np.add(adata_norm,random_matrix)

# new anndata
new_adata = sc.AnnData(new_matrix)#, dtype=X.dtype)
new_adata.var = adata.var # append metadata
new_adata.obs = adata.obs # append metadata


## UMAP and clustering for both the files

