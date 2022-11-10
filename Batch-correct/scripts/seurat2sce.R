#!/home/jacopo/miniconda3/envs/seurat_env/bin/Rscript

library(Seurat)

file <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]

dat <- readRDS(file)

dat <- as.SingleCellExperiment(dat)

saveRDS(dat, output)