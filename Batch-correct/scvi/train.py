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
parser.add_argument("-t","--test",
                    dest='test',
                    type=bool,
                    default=False,
                    help='Run test parameters. Default False')


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
import optuna
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
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 10,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1
}
## Running params

if args.test:                           
    scvi_epochs = 10
    scanvi_epochs = 5
    trials=10
    jobs=2
    batch_size=4096
    hvg=4000
    fract=.01
else:
    scvi_epochs = 1000
    scanvi_epochs = 100
    trials=50
    jobs=4
    batch_size=4096
    hvg=4000
    fract=.1



cpus=os.cpu_count()
n_cores=jobs//cpus

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
def objective(trial,adata_og,adata,batch_size,label_key,batch_key):
    # sample hyperparameters
    n_hidden = trial.suggest_int("n_hidden", 64, 512, step=64)
    n_latent=trial.suggest_int("n_latent", 10, 50, step=5)
    n_layers = trial.suggest_int("n_layers",1,5)
    lr=trial.suggest_float("lr", 1e-5, 1e-3, log=True)


    # setup scVI
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key,labels_key=label_key)
    model = scvi.model.SCVI(adata, 
                            n_hidden=n_hidden, 
                            n_latent=n_latent, 
                            n_layers=n_layers)
    model.train(max_epochs=scvi_epochs, 
                batch_size=batch_size,
                plan_kwargs={"lr":lr},
                early_stopping=early_stopping_kwargs,
                load_sparse_tensor=True,
                train_size=.8,
                validation_size=.1)
  

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        model,
        adata=adata,
        labels_key=label_key,
        unlabeled_category="Unknown",
    )


    scanvi_model.train(max_epochs=scanvi_epochs, 
                        batch_size=batch_size,
                        early_stopping=early_stopping_kwargs,
                        train_size=.8,
                        validation_size=.1)
    
    adata_int=compile_red_dim(adata,model=model)
    # get score
    score = evaluate_model(adata,adata_int,label_key,batch_key)
    trial.set_user_attr("model", model)

    return score



### define metrics
def evaluate_model(adata,adata_int,label_key,batch_key):
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
                            n_cores=n_cores)
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
                            n_cores=n_cores)
    silouhette_batch=scib.metrics.silhouette_batch(adata_int,
                                                    batch_key=batch_key,
                                                    label_key=label_key,
                                                    embed='X_umap',
                                                    verbose=False)
    ### combine biological scores
    bio_scores=np.mean([ari,clisi,isolated_labels_asw,nmi,silouhette])
    ### combine batch correction scores
    batch_scores=np.mean([graph_connectivity,ilisi,silouhette_batch])                
    return np.mean([ari,clisi,isolated_labels_asw,nmi,silouhette,graph_connectivity,ilisi,silouhette_batch])


if __name__ == "__main__":
    adata = sc.read(args.input_object) # read in data
    sc.pp.filter_cells(adata, min_counts=100) # filter cells
    sc.pp.subsample(adata, fraction=fract) # subsample data to fast parameter search
    sc.pp.filter_genes(adata, min_counts=1) # filter genes
    adata=preprocessing(adata, hvg) # preprocess data (find hvg)
    adata_og=compile_red_dim(adata) # compile reduced dimensions
    study = optuna.create_study(direction="maximize")
    study.optimize(lambda trial: objective(trial,
                                        adata_og,
                                        adata,
                                        batch_size,
                                        args.label_key,
                                        args.batch_key), 
                    n_trials=trials, 
                    n_jobs=jobs)

    best_trial = study.best_trial
    # best_model = best_trial.user_attrs["model"]
    # best_model.save(args.output_prefix)
    best_params = study.best_trial.params
    with open(f"{args.output_prefix}_params.json", "w") as f:
        json.dump(best_params, f)
    fig = optuna.visualization.plot_optimization_history(study)
    fig.write_image(f"{args.output_prefix}_optuna_opt.png")
    print(f"Best trial: {best_trial.number}")
    print(f"Best parameters: {best_trial.params}")
    print(f"Best score: {best_trial.value}")
    