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


## Running params

hvg=4000
cpus=os.cpu_count()


### define functions
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


def compile_red_dim(adata, model=None, resolution=None):
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
    if resolution == None:
        resolution=scib.metrics.cluster_optimal_resolution(adata_tmp, label_key=args.label_key, cluster_key='opt_cluser')[0]# get optimal resolution
    else:
        pass
    sc.tl.leiden(adata_tmp,resolution=resolution)
    sc.tl.umap(adata_tmp, min_dist=0.3, copy=False)
    return adata_tmp



def evaluate_model(adata_int,label_key,batch_key):
    """
    Evaluate the model
    Parameters
    ----------
    adata
        AnnData object
    adata_int
        AnnData object with integrated dimensions
    label_key
        Key in adata.obs that defines cell types
    batch_key
        Key in adata.obs that defines batches
    Returns
    -------
    float
        Sum of the metrics
    """
    # biological conservation:
    ari=scib.metrics.ari(adata_int,
                        cluster_key='leiden',
                        label_key=label_key)
    clisi=scib.metrics.clisi_graph(adata_int,
                            label_key=label_key, 
                            type_='embed', 
                            use_rep='X_umap',
                            n_cores=cpus)
    isolated_labels_asw=scib.metrics.isolated_labels_asw(adata_int, 
                                                        label_key=label_key, 
                                                        batch_key=batch_key,
                                                        embed='X_umap',
                                                        verbose=False)
    nmi=scib.metrics.nmi(adata_int,
                        cluster_key='leiden',
                        label_key=label_key)
    silouhette=scib.metrics.silhouette(adata_int,
                                        label_key=label_key,
                                        embed='X_umap')
    # batch correction
    graph_connectivity=scib.metrics.graph_connectivity(adata_int,
                                                        label_key=label_key)
    ilisi=scib.metrics.ilisi_graph(adata_int,
                            batch_key=batch_key,
                            type_='embed',
                            use_rep='X_umap',
                            n_cores=cpus)
    silouhette_batch=scib.metrics.silhouette_batch(adata_int,
                                                    batch_key=batch_key,
                                                    label_key=label_key,
                                                    embed='X_umap',
                                                    verbose=False)
    metrics={'ari':ari,
            'clisi':clisi,
            'isolated_labels_asw':isolated_labels_asw,
            'nmi':nmi,
            'silouhette':silouhette,
            'graph_connectivity':graph_connectivity,
            'ilisi':ilisi,
            'silouhette_batch':silouhette_batch}
    metrics={k:float(v) for k,v in metrics.items()}
    return metrics


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
    pred=model.posterior_predictive_sample(adata)
    pred=pred.to_scipy_sparse()
    mse=np.mean((adata.X-pred).power(2))
    rmse=np.sqrt(mse)
    mae=np.mean(np.abs(adata.X-pred))
    features_errors=np.mean(np.abs(adata.X-pred),axis=0)
    obs_errors=np.mean(np.abs(adata.X-pred),axis=1)
    errors={'mse':float(mse),
            'rmse':float(rmse),
            'mae':float(mae),
            'features_errors':features_errors.astype(float).tolist()[0],
            'obs_errors':obs_errors.flatten().astype(float).tolist()[0]}
    return errors


if __name__ == "__main__":
    adata = sc.read(args.input_object) # read in data
    adata=preprocessing(adata,hvg) # preprocess data
    model=scvi.model.SCVI.load('results/BM_dataset_scvi_model/',adata) # initialize model
    adata_int=compile_red_dim(adata,model=model)
    scvi_metrics = evaluate_model(adata_int,args.label_key,args.batch_key)
    scvi_errors=reconstruction_error(model,adata) # get reconstruction error
    with open(f"{args.output_prefix}_scvi_metrics.json", "w") as f:
        json.dump(scvi_metrics, f) # save
    with open(f"{args.output_prefix}_scvi_errors.json", "w") as f:
        json.dump(scvi_errors, f)
    adata_int.write(f"{output_prefix}_scvi_corrected.h5ad") # save corrected data
    scanvi_model=scvi.model.SCANVI.load('results/BM_dataset_scanvi_model/',adata) # initialize model
    adata_scanvi=compile_red_dim(adata,model=scanvi_model)
    scanvi_metrics = evaluate_model(adata_scanvi,args.label_key,args.batch_key)
    scanvi_errors=reconstruction_error(scanvi_model,adata) # get reconstruction error
    # save scores and results as json

    with open(f"{args.output_prefix}_scanvi_metrics.json", "w") as f:
        json.dump(scanvi_metrics, f)
    with open(f"{args.output_prefix}_scanvi_errors.json", "w") as f:
        json.dump(scanvi_errors, f)
    adata_scanvi.write(f"{output_prefix}_scanvi_corrected.h5ad") # save corrected data

    adata_og=compile_red_dim(adata)
    og_metrics = evaluate_model(adata_og,args.label_key,args.batch_key)
    with open(f"{args.output_prefix}_og_metrics.json", "w") as f:
        json.dump(og_metrics, f)
    adata_og.write(f"{output_prefix}_og.h5ad") # save corrected data

