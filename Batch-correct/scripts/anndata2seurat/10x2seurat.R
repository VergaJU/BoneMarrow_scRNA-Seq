#!/home/jacopo/miniconda3/envs/seurat_env/bin/Rscript

library(Seurat)
library(SingleCellExperiment)

dat <- Read10X("temp/")
metadata <- read.csv("./temp/metadata.csv", row.names="X")
dat <- CreateSeuratObject(dat, meta.data=metadata)
umap <- read.csv("./temp/umap.csv", row.names="X")
dat <- as.SingleCellExperiment(dat)
reducedDim(dat, "umap") <- umap
dat <- as.Seurat(dat)

saveRDS(dat, "nk_raw_umap_red.Rds")