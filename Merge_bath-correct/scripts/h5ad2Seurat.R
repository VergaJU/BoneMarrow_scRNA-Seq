#!/home/jacopo/anaconda3/envs/convert/bin/Rscript
library(reticulate)
library(Seurat)
library(Matrix)
library(optparse)

option_list = list(
  make_option(
    c("-i", "--input_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to input file.'
  ),
  make_option(
    c("-o", "--output_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to output file.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))


#make sure you have python installed. This works with Python 3.8.5
#if you don't already have it, install the python "scanpy" package using reticulate use the following to install it:
#py_install("scanpy")
## method find here https://github.com/satijalab/seurat/issues/4711

sc <- import("scanpy")

adata <- sc$read_h5ad(opt$input_file)

counts <- t(adata$X)
colnames(counts) <-  adata$obs_names$to_list()
rownames(counts) <-  adata$var_names$to_list()
counts <- Matrix::Matrix(as.matrix(counts), sparse = T)
#seurat <- as(counts, "dgCMatrix")

seurat <- CreateSeuratObject(counts)
seurat <- AddMetaData(seurat,  adata$obs)

saveRDS(seurat, file = opt$output_file)