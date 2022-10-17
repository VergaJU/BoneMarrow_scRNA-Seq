#!/home/jacopo/miniconda3/envs/harmony/bin/Rscript

library(harmony)
library(Seurat)
library(optparse)
library(stringr)


option_list = list(
  make_option(
    c("-i", "--input_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to seurat object file, already nurmalized and PCA run.'
  ),
  make_option(
    c("-b", "--batch_key"),
    action = "store",
    default = "batch",
    type = 'character',
    help = 'Variable where to store the batch.'
  ),
  make_option(
    c("-c", "--celltype_key"),
    action = "store",
    default = "label",
    type = 'character',
    help = 'Variable where to store the batch.'
  ),
  make_option(
    c("-t", "--theta_batch"),
    action = "store",
    default = 2,
    type = 'double',
    help = 'Diversity clustering penalty parameter. Higher theta more diverse clusters, theta=0 no diversity. Default = 2'
  ),
  make_option(
    c("-l", "--theta_celltype"),
    action = "store",
    default = 0,
    type = 'double',
    help = 'Diversity clustering penalty parameter. Higher theta more diverse clusters, theta=0 no diversity. Default = 0'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suffix <- paste0("_harmony_",opt$theta_batch,"_",opt$theta_celltype,".Rds") 

output <-str_replace(opt$input_file, ".Rds", suffix)

dat <- readRDS(opt$input_file)

vars = c(opt$batch_key, opt$celltype_key)
thetas = c(opt$theta_batch, opt$theta_celltype)

dat <- RunHarmony(dat, vars, theta = thetas, max.iter.harmony = 100)
saveRDS(dat, output)