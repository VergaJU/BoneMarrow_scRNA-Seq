#!/home/jacopo/miniconda3/envs/seurat_env/bin/Rscript

# TODO: improve docstring
library(Seurat)
library(SingleCellExperiment)
library(stringr)

file_list <- commandArgs(trailingOnly = TRUE)[1]
condition_name <- commandArgs(trailingOnly = TRUE)[2]

merge.data <- function(filepath,condition){
    # Load file paths
    files <- read.csv(filepath, sep = "\t", header=F)
    # create variables
    dats = list()
    labels = list()
    samples = vector()

    # append Rds, labels and sample names
    for(i in 1:nrow(files)){
        rds = files[i,1]
        labs = files[i,2]
        dats[i] = readRDS(rds) # Load Rds
        labels[[i]] = read.csv(file=(labs),sep="\t") # Load celltype labels
        samples[i] = str_extract(str_extract(files[i,1], "[^/]+$"), "[A-Z]+[0-9]+") # load sample name
        dats[[i]] = AddMetaData(dats[[i]],labels[[i]])
        dats[[i]]$orig.ident = samples[i]
    }

    # merge datas
    merged.dat = merge(dats[[1]], y = dats[2:length(dats)], add.cell.ids = samples)
    
    # Save seurat Object
    saveRDS(merged.dat, file = paste0(condition, "_dataset_seurat.Rds"))
    merged.dat = as.SingleCellExperiment(merged.dat)
    saveRDS(merged.dat, file = paste0(condition, "_dataset_sce.Rds"))

}

merge.data(file_list, condition_name)