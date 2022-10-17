#!/home/jacopo/miniconda3/envs/seurat_env/bin/Rscript

library(Seurat)
library(Matrix)

args <- commandArgs(trailingOnly = T) # save argument in a vector

#This set of functions exports the counts of Seurat to matrixmarket
#Note that only Gene Expression is exported

Seurat_to_MM10X <- function(Seurat,prefix){
    writeMM(Seurat@assays$RNA@counts,paste(prefix,"matrix.mtx",sep="")) # Write matrix w/ counts
    write(x = rownames(Seurat@assays$RNA@counts),file = paste(prefix,"genes.tsv",sep="")) # write gene
    ## TODO: genes should have 2 columns, TO BE FIXED
    write(x = colnames(Seurat@assays$RNA@counts),file = paste(prefix,"barcodes.tsv",sep="")) # write barcodes
    write.table(x = Seurat@meta.data, file = paste(prefix,"metadata.tsv",sep=""), # write metadata
    quote=FALSE,sep="\t")
}

dat <- readRDS(args[1]) # Read input file

Seurat_to_MM10X(dat,"temp/") # convert the file in the "temp" directory