#!/home/jacopo/miniconda3/envs/seurat_env/bin/Rscript

library(Seurat)
library(SingleCellExperiment)
library(Matrix)

dat <- Read10X("temp/")
metadata <- read.csv("./temp/metadata.csv", row.names="X")
dat <- CreateSeuratObject(dat, meta.data=metadata)
umap <- read.csv("./temp/umap.csv", row.names="X")
pca <- read.csv("./temp/pca.csv", row.names = "X")
dat <- as.SingleCellExperiment(dat)
reducedDim(dat, "umap") <- umap
reducedDim(dat, "pca") <- pca
dat <- as.Seurat(dat)
hvg <- read.csv("./temp/hvg.csv", header = F)
hvg <- unlist(c(hvg[1]))

saveRDS(dat, "DE_data.Rds")

