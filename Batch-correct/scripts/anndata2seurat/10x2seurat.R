#!/home/jacopo/miniconda3/envs/seurat_env/bin/Rscript


suppressMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(Matrix)
}
)

args <- commandArgs(trailingOnly = T) # save argument in a vector

dat <- Read10X("temp/")
metadata <- read.csv("./temp/metadata.csv", row.names="X")
dat <- CreateSeuratObject(dat, meta.data=metadata)


saveRDS(dat,args[1])

