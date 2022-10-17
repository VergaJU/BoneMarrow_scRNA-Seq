#!/home/jacopo/anaconda3/envs/seurat_env/bin/Rscript

library(Seurat)

file <- commandArgs(trailingOnly = TRUE)[1]

dat <- readRDS(file)

dat <- as.SingleCellExperiment(dat)

saveRDS(dat, file)