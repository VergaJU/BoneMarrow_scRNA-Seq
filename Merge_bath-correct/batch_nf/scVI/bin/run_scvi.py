#!/opt/conda/bin/python3


import argparse
import scanpy as sc
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata
import scipy

# read in parameters
parser = argparse.ArgumentParser(description='Input/Output files and method parameters')
parser.add_argument("--input_object",
                    dest='input_object',
                    type=str,
                    help ='Input h5ad file.')
parser.add_argument("--output_prefix",
                    dest='output_prefix',
                    type=str,
                    help='Name to use as prefix for saving objects.')
parser.add_argument("--batch_key",
                    dest='batch_key',
                    type=str,
                    default='orig.ident',
                    help ='Obs key defining batch.')
parser.add_argument("--celltype_key",
                    dest="celltype_key",
                    type=str,
                    default="label",
                    help = "Obs key defining cell_types")
parser.add_argument('--hvg',
                    dest='hvg',
                    type=int,
                    help='Number of HVGs.')
parser.add_argument('--latent',
                    dest='latent',
                    type=int,
                    help='Number of latent dimensions')

args = parser.parse_args()

    
# setup fixed parameters
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)
early_stopping_kwargs = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 10,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1
}

# read in data 
adata = sc.read(args.input_object)

# convert the sparse matrix from csc to csr to speed up the training
matrix = adata.X
matrix = scipy.sparse.csr_matrix(matrix)
adata_csr = anndata.AnnData(matrix_csr)
adata_csr.obs = adata.obs
adata_csr.var = adata.var
adata_csr.obsm = adata.obsm
adata = adata_csr

# subset to variable genes
sc.pp.log1p(adata) # copute log
sc.pp.highly_variable_genes(adata, n_top_genes=args.hvg, batch_key=args.batch_key, subset=True)

# get counts
#adata.X = adata[:, adata.var_names].layers['counts']

# setup scarches SCVI
sca.models.SCVI.setup_anndata(adata, batch_key=args.batch_key)
vae = sca.models.SCVI(adata,
        n_layers=2,
        n_latent=args.latent,
        encode_covariates=True,
        deeply_inject_covariates=False,
        use_layer_norm="both",
        use_batch_norm="none"
        )

# train model
print('Training model...')
vae.train(max_epochs=500, early_stopping=early_stopping_kwargs)
print('Done!')

# save model
ref_path = 'ref_model_scVI' + args.output_prefix
vae.save(ref_path, overwrite=True)
reference_latent = sc.AnnData(vae.get_latent_representation())
reference_latent.obs_names = adata.obs_names
reference_latent.obs = adata.obs
save_path = args.output_prefix + '_latent.h5ad'
reference_latent.write(save_path)
