#! /home/jacopo/miniconda3/envs/scanpy_env/bin/python
import scipy.sparse as sparse
import scipy.io as sio
import scanpy as sc
import pandas as pd
import sys
import os
filename = sys.argv[1]

adata = sc.read(filename)

mtx = adata.X.transpose()

if not os.path.exists("./temp"):
    os.makedirs("./temp")

sio.mmwrite("./temp/matrix.mtx",mtx)
barcodes = pd.DataFrame(adata.obs_names)
barcodes.to_csv("./temp/barcodes.tsv", sep="\t", index=False, header=False)
genes= pd.DataFrame(adata.var_names)
genes[1] = genes[0]
genes.to_csv("./temp/genes.tsv", sep="\t", index=False, header=False)
adata.obs.to_csv("./temp/metadata.csv")
