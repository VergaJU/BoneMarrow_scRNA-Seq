#! /home/jacopo/miniconda3/envs/scanpy_env/bin/python
import scipy.sparse as sparse
import scipy.io as sio
import scanpy as sc
import pandas as pd


adata = sc.read("MM_atlas_umap.h5ad")

mtx = adata.X.transpose()
sio.mmwrite("./temp/matrix.mtx",mtx)
barcodes = pd.DataFrame(adata.obs_names)
barcodes.to_csv("./temp/barcodes.tsv", sep="\t", index=False, header=False)
genes= pd.DataFrame(adata.var_names)
genes[1] = genes[0]
genes.to_csv("./temp/genes.tsv", sep="\t", index=False, header=False)
adata.obs.to_csv("./temp/metadata.csv")
pd.DataFrame(adata.obsm["X_umap"], columns = ["UMAP_1", "UMAP_2"], index=adata.obs_names).to_csv("./temp/umap.csv")
pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names).to_csv("./temp/pca.csv")
#hvg = adata.var.highly_variable
#hvg = pd.Series(hvg[hvg == True].index)
#hvg.to_csv("./temp/hvg.csv", index=None, header=None)
