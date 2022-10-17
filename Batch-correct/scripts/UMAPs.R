#!/home/jacopo/miniconda3/envs/seurat_env/bin/Rscript

library(Seurat)
library(stringr)
library(tidyverse)

file <- commandArgs(trailingOnly = TRUE)[1] # Seurat File
condition <- commandArgs(trailingOnly = TRUE)[2] # Condition name
batch <- commandArgs(trailingOnly = TRUE)[3] # Batch (yes or no)

dat <- readRDS(file) # load rds file
dat <- NormalizeData(dat) # normalize counts
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 4000) # find var features to compute neighbors
#all.genes <- rownames(dat) # list of genes
dat <- ScaleData(dat) # scale counts (var 1, avg 0)
dat <- RunPCA(dat, features = VariableFeatures(object = dat), verbose = F, seed.use = 1) # PCA using variable features
dat <- FindNeighbors(dat, dims = 1:15) # compute neighbors tree with first 15 PCs
dat <- RunUMAP(dat, dims = 1:15, seed.use = 1) # compute UMAP with first 15 PCs

svg(paste0(condition,"_celltypes_", batch, "-batch.svg")) # open svg object
# plot UMAP, colors grouped by label (celltype)
DimPlot(dat, reduction = 'umap', group.by="label", label=T, label.size=3, raster=FALSE, seed = 1)
# close device, save svg
dev.off()


svg(paste0(condition,"_samples_",batch, "-batch.svg")) # open svg object
# plot UMAP, colors grouped by orig.ident (sample)
DimPlot(dat, reduction = 'umap', group.by="orig.ident", label=T, label.size=3, raster=FALSE, seed = 1)
# close device, save svg
dev.off()