#!/usr/bin/python3

import argparse 



# read in parameters
parser = argparse.ArgumentParser(description='Input/Output files and method parameters')
parser.add_argument("-i", "--input_object",
                    dest='input_object',
                    type=str,
                    help ='Input h5ad file.')
parser.add_argument("-o", "--output_prefix",
                    dest='output_prefix',
                    type=str,
                    help='Name to use as prefix for saving objects.')
parser.add_argument("-b", "--batch_key",
                    dest='batch_key',
                    type=str,
                    default='batch',
                    help ='Obs key defining batch.Default batch')
parser.add_argument("-l", "--label",
                    dest="label_key",
                    type=str,
                    default='label',
                    help='Label key. Default label')
parser.add_argument("-p", "--params",
                    dest="params",
                    type=str,
                    help='JSON file with parameters')


args = parser.parse_args()


import os
import tempfile

import scanpy as sc
import scvi
import seaborn as sns
import torch
from rich import print
import scib
import numpy as np
from plotly.io import show
import json
import warnings
warnings.filterwarnings("ignore")

# Set parameters for figures and torch
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()

torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)
early_stopping_kwargs = {
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_min_delta": 0
}
## Running params

scvi_epochs = 1000
scanvi_epochs = 100
batch_size=16384
hvg=4000



cpus=os.cpu_count()

def get_params(params):
    with open(params) as f:
        data = json.load(f)
    return data

### define functions
# preprocessing
## adata_pp is processed
def preprocessing(adata, high_variable_genes)->sc.AnnData:
    """
    Preprocesses the AnnData object for scVI training.
    Parameters
    ----------
    adata
        AnnData object
    high_variable_genes
        Number of high variable genes to select
    batch_key
        Key in adata.obs that defines batches

    Returns
    -------
    AnnData object
        adata with preprocessed counts
    """
    # Saving count 
    adata.raw=adata
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=high_variable_genes,
        subset=True,
    )
    return adata


# compile reduced dimensions
def compile_red_dim(adata, model=None):
    """
    Compiles the reduced dimensions of the AnnData object.
    Parameters
    ----------
    adata
        AnnData object
    model_key
        Key to store the model
    model
        Model to use for the reduced dimensions
    Returns
    -------
    AnnData object
        adata with reduced dimensions
    """
    adata_tmp=adata.copy()
    if model == None:
        model_key='pca'
        sc.pp.normalize_total(adata_tmp)
        sc.pp.log1p(adata_tmp)
        sc.pp.pca(adata_tmp)
    else:
        adata_tmp.obsm['X_pca'] = model.get_latent_representation() # if latendt dimensions from a model, add them as X_pca
    sc.pp.neighbors(adata_tmp, use_rep='X_pca')
    resolution=scib.metrics.cluster_optimal_resolution(adata_tmp, label_key=args.label_key, cluster_key='opt_cluser')[0]# get optimal resolution
    # resolution=.6
    sc.tl.leiden(adata_tmp,resolution=resolution)
    sc.tl.umap(adata_tmp, min_dist=0.3, copy=False)
    return adata_tmp



## objective function for Optuna
def train(adata,batch_size,label_key,batch_key,params):
    data = get_params(params)
    # sample hyperparameters
    n_hidden = data['n_hidden']
    n_latent=data['n_latent']
    n_layers = data['n_layers']
    lr=data['lr']


    # setup scVI
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key,labels_key=label_key)
    model = scvi.model.SCVI(adata, 
                            n_hidden=n_hidden, 
                            n_latent=n_latent, 
                            n_layers=n_layers)
    model.train(max_epochs=scvi_epochs, 
                batch_size=batch_size,
                plan_kwargs={"lr":lr},
                early_stopping=True,
                load_sparse_tensor=True,
                train_size=.8,
                validation_size=.1,
                **early_stopping_kwargs)
      

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        model,
        adata=adata,
        labels_key=label_key,
        unlabeled_category="Unknown",
    )


    scanvi_model.train(max_epochs=scanvi_epochs, 
                        batch_size=batch_size,
                        plan_kwargs={"lr":lr/10},
                        train_size=.8,
                        validation_size=.1)
    


    return model,scanvi_model






def reconstruction_error(model, adata):
    """
    Calculate the reconstruction error
    Parameters
    ----------
    model
        scVI model
    adata
        AnnData object
    Returns
    -------
    dict
        Dictionary with the reconstruction error
    """
    pred=model.posterior_predictive_sample(adata.X)
    mse=np.mean((adata.X-pred).powzer(2))
    rmse=np.sqrt(mse)
    mae=np.mean(np.abs(adata.X-pred))
    features_errors=np.mean(np.abs(adata.X-pred),axis=0)
    obs_errors=np.mean(np.abs(adata.X-pred),axis=1)
    errors={'mse':mse,
            'rmse':rmse,
            'mae':mae,
            'features_errors':features_errors,
            'obs_errors':obs_errors}
    return errors



if __name__ == "__main__":
    adata = sc.read(args.input_object) # read in data
    print(f"Reading data from {args.input_object}")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    adata=preprocessing(adata, hvg) # preprocess data (find hvg)
    print(f"Hvg selected: {adata.n_vars}")
    print(f"Training scVI model with {scvi_epochs} epochs")
    model,scanvi_model=train(adata,batch_size,args.label_key,args.batch_key,args.params) # train model

    print("Training complete, saving objects")
    model.save(f"{args.output_prefix}_scvi_model") # save model
    scanvi_model.save(f"{args.output_prefix}_scanvi_model") # save model
