#!/home/jacopo/miniconda3/envs/swarms/bin/python3
import pandas as pd
import numpy as np
import scanpy as sc
import sys
#import pyswarms as ps
#import scib

inputfile = sys.argv[1]
outputfile = inputfile[:-5] + "_clustered.h5ad"
outputfile2 = inputfile[:-5] + "_new_clustered.h5ad"

# Read file
adata = sc.read(inputfile)


# filter genes expressed in less than 3 cells
sc.pp.filter_genes(adata, min_cells = 3)


# normalize data
sc.pp.normalize_total(adata, target_sum=1e4)

# Run UMAP and clustering as reference
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)
adata.write(outputfile)

# new matrix used to be modified by birds
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


## UMAP and clustering of the new anndata

sc.pp.log1p(new_adata)
sc.pp.highly_variable_genes(new_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(new_adata, svd_solver='arpack')
sc.pp.neighbors(new_adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(new_adata)
sc.tl.leiden(new_adata)
new_adata.write(outputfile)