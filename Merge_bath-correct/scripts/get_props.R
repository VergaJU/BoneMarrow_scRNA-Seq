#!/home/jacopo/miniconda3/envs/seurat_env/bin/Rscript

library(Seurat)
library(tidyverse)
library(stringr)

dataset <- commandArgs(trailingOnly = TRUE)[1] # input: path/to/seurat/object

get_proportions = function(dataset){
    dat = readRDS(dataset) # load seurat object
    cells = table(dat$label) # store cell frequencies in a new variablr
    cells = as.data.frame(cells) # convert variable in dataframe
    colnames(cells) = c("celltype","freq") # rename columns in a meaningful way
    total = sum(cells$freq) # get total number of cells to compute proportions
    cells = cells %>%   
        mutate(proportions = freq/total) %>% # get proportion of celltype
        mutate(perc = proportions*100) # get percentage

    write_csv(cells, paste0(str_extract(dataset, "[^.]+"), "_proportions.csv")) # save file
}

get_proportions(dataset)