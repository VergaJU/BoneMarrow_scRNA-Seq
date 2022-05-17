#!/home/jacopo/miniconda3/envs/seurat_env/bin/Rscript
library(Seurat)
library(stringr)
library(tidyverse)

file <- commandArgs(trailingOnly = TRUE)[1] # Seurat File
condition <- commandArgs(trailingOnly = TRUE)[2] # Condition name
batch <- commandArgs(trailingOnly = TRUE)[3] # Batch (yes or no)

dat <- readRDS("./SMM_dataset_seurat.Rds")
dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(dat)
dat <- ScaleData(dat)
dat <- RunPCA(dat, features = VariableFeatures(object = dat), verbose = F, seed.use = 1)
dat <- FindNeighbors(dat, dims = 1:15)
dat <- RunUMAP(dat, dims = 1:15, seed.use = 1)

svg(paste0(condition,"_celltypes_", batch, "-batch.svg"))
DimPlot(dat, reduction = 'umap', group.by="label", label=T, label.size=3, seed = 1)
dev.off()


svg(paste0(condition,"_samples_",batch, "-batch.svg"))
DimPlot(dat, reduction = 'umap', group.by="orig.ident", label=T, label.size=3, seed = 1)
dev.off()