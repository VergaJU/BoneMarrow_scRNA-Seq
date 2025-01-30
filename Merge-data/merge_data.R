#!/home/jacopo/anaconda3/envs/seurat_env/bin/Rscript


library(optparse)




option_list = list(
  make_option(
    c("-i", "--input_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to input file, it must be a file containing the paths for the files to be merged and metadata.'
  ),
  make_option(
    c("-o", "--output_prefix"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Prefix for naming output file.'
  ),
  make_option(
    c("-b", "--batch_column"),
    action = "store",
    default = "batch",
    type = 'character',
    help = 'Column with batch names..'
  ),
  make_option(
    c("-p", "--path_column"),
    action = "store",
    default = "file",
    type = 'character',
    help = 'Column with the paths for the files '
  )
)

opt <- parse_args(OptionParser(option_list=option_list))


library(Seurat)
library(SingleCellExperiment)
library(stringr)
library(parallel)
library(BiocParallel)

multicoreParam <- MulticoreParam(workers = detectCores())

# merge data inputted
## input list of files and condition
## output: seurat and sce objects named as condition
merge.data <- function(input_file,batch_column,path_column,output_prefix){
    # Load file paths
    df <- read.csv(input_file, header=T)
    files <- df[path_column]
    batches <- df[batch_column]
    meta <- setdiff(colnames(df), c(batch_column,path_column))
    # create variables
    dats = list()
    cat("Merging", length(files[,1]), "files\nBatch key:", batch_column, "\n")

    # append Rds, labels and sample names
    for(i in 1:nrow(files)){
        rds = files[i,1] # load path for rds file
        batch=batches[i,1] # get batch
        dats[i] = readRDS(rds) # Load Rds
        dats[[i]]$batch = batch # add batch metadata
        for(m in meta){
          dats[[i]]@meta.data[m] <- df[i,m] # add rest of metadata
        }
    }

    # merge datas
    # merge(first obj, y = <vector with other objects>, add.cell.ids = <vector of samples names>)
    merged.dat = merge(dats[[1]], y = dats[2:length(dats)], add.cell.ids = batches$batch, BPPARAM=multicoreParam)
    
    # Save seurat Object
    saveRDS(merged.dat, file = paste0(output_prefix, "_dataset_seurat.Rds"))
    # convert to sce object
    merged.dat = as.SingleCellExperiment(merged.dat)
    # fix assay name and remove logcounts
    #assay(merged.dat, "RNA") <- assay(merged.dat, "counts")
    #assay(merged.dat, "counts") <- NULL
    #assay(merged.dat, "logcounts") <- NULL
    # save sce object
    saveRDS(merged.dat, file = paste0(output_prefix, "_dataset_sce.Rds"))

}


merge.data(opt$input_file, opt$batch_column, opt$path_column, opt$output_prefix)
