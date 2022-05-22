#!/home/jacopo/miniconda3/envs/seurat_env/bin/Rscript


library(Seurat)
library(SingleCellExperiment)
library(stringr)

 # tsv file containing paths to the files to be merged and cell labels obtained with stringR
 # path/to/file path/to/labels
file_list <- commandArgs(trailingOnly = TRUE)[1]
# condition to rename the merged file (string)
condition_name <- commandArgs(trailingOnly = TRUE)[2] 


# merge data inputted
## input list of files and condition
## output: seurat and sce objects named as condition
merge.data <- function(filepath,condition){
    # Load file paths
    files <- read.csv(filepath, sep = "\t", header=F)
    # create variables
    dats = list()
    labels = list()
    samples = vector()

    # append Rds, labels and sample names
    for(i in 1:nrow(files)){
        rds = files[i,1] # load path for rds file
        labs = files[i,2] #  load path to labels
        dats[i] = readRDS(rds) # Load Rds
        labels[[i]] = read.csv(file=(labs),sep="\t") # Load celltype labels
        samples[i] = str_extract(str_extract(files[i,1], "[^/]+$"), "[A-Z]+[0-9]+") # load sample name
        dats[[i]] = AddMetaData(dats[[i]],labels[[i]]) # add cell labels to seurat object
        dats[[i]]$orig.ident = samples[i] # rename orig ident with sample name
    }

    # merge datas
    # merge(first obj, y = <vector with other objects>, add.cell.ids = <vector of samples names>)
    merged.dat = merge(dats[[1]], y = dats[2:length(dats)], add.cell.ids = samples)
    
    # Save seurat Object
    saveRDS(merged.dat, file = paste0(condition, "_dataset_seurat.Rds"))
    # convert to sce object
    merged.dat = as.SingleCellExperiment(merged.dat)
    # fix assay name and remove logcounts
    assay(dat, "RNA") <- assay(dat, "counts")
    assay(dat, "counts") <- NULL
    assay(dat, "logcounts") <- NULL
    # save sce object
    saveRDS(merged.dat, file = paste0(condition, "_dataset_sce.Rds"))

}


merge.data(file_list, condition_name)