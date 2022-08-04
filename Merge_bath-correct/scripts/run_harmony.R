#!/home/jacopo/anaconda3/envs/harmony/bin/Rscript

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
    c("-t", "--theta"),
    action = "store",
    default = 2,
    type = 'double',
    help = 'Diversity clustering penalty parameter. Higher theta more diverse clusters, theta=0 no diversity. Default = 2'
  ),
  make_option(
    c("-l", "--lambda"),
    action = "store",
    default = 1,
    type = 'double',
    help = 'Ridge regression penalty parameter. Must be positive, smaller more aggressive correction. Default = 1'
  ),
  make_option(
    c("-s", "--sigma"),
    action = "store",
    default = 0.1,
    type = 'double',
    help = 'Width of soft kmeans clusters. Smaller values more similar to hard clustering. Default = 0.1'
  ),
  make_option(
    c("-b", "--batch_key"),
    action = "store",
    default = "batch",
    type = 'character',
    help = 'Variable where to store the batch.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suffix <- paste0("_harmony_",opt$theta,"_",opt$lambda,"_",opt$sigma,".Rds") 

output <-str_replace(opt$input_file, ".Rds", suffix)

dat <- readRDS(opt$input_file)

dat <- RunHarmony(dat, opt$batch_key, theta = opt$theta, lambda = opt$lambda, sigma = opt$sigma, max.iter.harmony = 100)
saveRDS(dat, output)