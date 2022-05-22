#!/home/jacopo/miniconda3/envs/seurat_env/bin/Rscript


library(Seurat)
library(SingleCellExperiment)
library(stringr)

 # tsv file containing paths to the files to be merged and cell labels obtained with stringR
 # path/to/file path/to/labels
file_list <- commandArgs(trailingOnly = TRUE)[1]
# condition to rename the merged file (string)
condition_name <- commandArgs(trailingOnly = TRUE)[2] 

# select 5 random files
files <- read.csv(file_list, sep = "\t", header=F)
x <- sample(1:length(files$V1), 10, replace=F)

# merge data inputted
## input list of files and condition
## output: seurat and sce objects named as condition
merge.data <- function(filepath,condition){
    # Load file paths
    # create variables
    dats = list()
    labels = list()
    samples = vector()

    # append Rds, labels and sample names
for(i in 1:length(x)){
    rds = files[x[i],1] # load path for rds file
    labs = files[x[i],2] #  load path to labels
    dats[i] = readRDS(rds) # Load Rds
    labels[[i]] = read.csv(file=(labs),sep="\t") # Load celltype labels
    samples[i] = str_extract(str_extract(files[x[i],1], "[^/]+$"), "[A-Z]+[0-9]+") # load sample name
    dats[[i]] = AddMetaData(dats[[i]],labels[[i]]) # add cell labels to seurat object
    dats[[i]]$orig.ident = samples[i] # rename orig ident with sample name
    dats[[i]] = subset(x=dats[[i]], downsample=500)
}

    # merge datas
    # merge(first obj, y = <vector with other objects>, add.cell.ids = <vector of samples names>)
    merged.dat = merge(dats[[1]], y = dats[2:length(dats)], add.cell.ids = samples)
    
    # Save seurat Object
    saveRDS(merged.dat, file = paste0(condition, "_dataset_seurat.Rds"))
    # convert to sce object
    merged.dat = as.SingleCellExperiment(merged.dat)
    # fix assay name and remove logcounts
    assay(merged.dat, "RNA") <- assay(merged.dat, "counts")
    assay(merged.dat, "counts") <- NULL
    assay(merged.dat, "logcounts") <- NULL
    # save sce object
    saveRDS(merged.dat, file = paste0(condition, "_dataset_sce.Rds"))

}


merge.data(file_list, condition_name)