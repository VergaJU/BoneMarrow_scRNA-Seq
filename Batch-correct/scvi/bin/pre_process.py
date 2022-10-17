#!/opt/conda/bin/python3


import argparse
import scanpy as sc
import pandas as pd
import scipy
import anndata

# read in parameters
parser = argparse.ArgumentParser(description='Input/Output files and method parameters')
parser.add_argument("--input_object",
                    dest='input_object',
                    type=str,
                    help ='Input h5ad file.')
parser.add_argument("--output_object",
                    dest='output_object',
                    type=str,
                    help='Name to use for preprocessed file.')
parser.add_argument("--batch_key",
                    dest='batch_key',
                    type=str,
                    default='orig.ident',
                    help ='Obs key defining batch.')
parser.add_argument('--hvg',
                    dest='hvg',
                    type=int,
                    help='Number of HVGs.')

args = parser.parse_args()

# read in data 
adata = sc.read(args.input_object)

# convert the sparse matrix from csc to csr to speed up the training
matrix = adata.X
matrix = scipy.sparse.csr_matrix(matrix)
adata_csr = anndata.AnnData(matrix)
adata_csr.obs = adata.obs
adata_csr.var = adata.var
adata_csr.obsm = adata.obsm
adata = adata_csr

# subset to variable genes
adata_log = adata.copy()
sc.pp.log1p(adata_log) # compute log
sc.pp.highly_variable_genes(adata_log, n_top_genes=args.hvg, batch_key=args.batch_key)
adata.var = adata_log.var
adata.uns["hvg"] = adata_log.uns["hvg"]
adata = adata[:, adata.var.highly_variable]

adata.write(args.output_object)
