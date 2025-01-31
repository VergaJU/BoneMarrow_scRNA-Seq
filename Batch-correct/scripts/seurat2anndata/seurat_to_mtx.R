#!/home/jacopo/miniconda3/envs/seurat_env/bin/Rscript

suppressMessages({
    library(Seurat)
    library(Matrix)
}
)

args <- commandArgs(trailingOnly = T) # save argument in a vector

#This set of functions exports the counts of Seurat to matrixmarket
#Note that only Gene Expression is exported

Seurat_to_MM10X <- function(dat,prefix){
    writeMM(dat@assays$RNA@counts,paste0(prefix,"matrix.mtx")) # Write matrix w/ counts
    gene_names <<- data.frame(V1=rownames(dat),V2=rownames(dat)) # get dataframe of gene names
    write.table(gene_names, paste0(prefix,"genes.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) # write genes
    write(x = colnames(dat@assays$RNA@counts),file = paste0(prefix,"barcodes.tsv")) # write barcodes
    write.table(x = dat@meta.data, file = paste0(prefix,"metadata.tsv"), quote=FALSE,sep="\t") # write metadata
}

dat <- readRDS(args[1]) # Read input file

print(dat)
Seurat_to_MM10X(dat,"temp/") # convert the file in the "temp" directory