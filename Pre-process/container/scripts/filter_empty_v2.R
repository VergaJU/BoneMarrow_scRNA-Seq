#!/opt/conda/envs/pre_process/bin/Rscript

library(tidyverse, quietly=T) # packages for data wrangling, visualization etc
library(DropletUtils, quietly=T)
library(scDblFinder, quietly=T)
library(Matrix, quietly=T)
library(tools, quietly = T)
library(R.utils, quietly = T)

set.seed(1) # set seed to increase reproducibility


# save arguments in a vector
args <- commandArgs(trailingOnly = T) 

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied\n 
       First argument: `path/to/files/`\n 
       Second argument (optional): FDR threshold ", call.=FALSE)
} else if (length(args)==1) {
  # default select FDR threshold of 0.1 file
  print("No FDR threshold supplied, default threshold of 0.1 selected")
  args[2] <- 0.1
}

filepath <- args[1]
files <- list.files(filepath)

if (file.exists(paste0(filepath, "barcodes.tsv.gz"))) {
  cat("\nFiles compressed, gunzipping...")
  for (f in files) {
    gunzip(paste0(filepath, f), remove = T)
    cat("...")
  }
  files <- list.files(filepath)
  cat("\nFiles decompressed, checking file names...")
} else {
  cat("\nFiles not compressed, checking file names...")
}


if ("genes.tsv" %in% files){
  genes.file <- "genes.tsv"
} else {
  genes.file <- "features.tsv"
}

cat("\nFilenames checked\nStart filtering empty drops with a threshold of", args[2],
    "\nThe process can take some minutes\n")

m <- Matrix::readMM(paste0(filepath, "/matrix.mtx")) # load matrix
genes  <- read.csv(paste0(filepath, genes.file), sep = "\t", header = F)# attach genes
rownames(m) <- genes[,2]
colnames(m) <- read.csv(paste0(filepath, "/barcodes.tsv"), sep = "\t", header = F)[,1] # attach barcodes
empty <- emptyDrops(m) # compute the probability of empty droplets
full <- empty$FDR <= args[2] # filter by selected FDR threshold
full[is.na(full)] <- FALSE # remove NAs
filtered_m <- m[,full]  # subset original matrix

cat("Empty droplets filtered correctly \n",
    "Starting filtering of doublets using scDblFinder\n")

dat <- SingleCellExperiment(assay = list(counts = filtered_m))
dat <- scDblFinder(dat)

singlets <- as.logical(as.integer(dat$scDblFinder.class) == 1)
filtered_m1 <- filtered_m[,singlets]  # subset original matrix
dim(m)
dim(filtered_m)
dim(filtered_m1)


DropletUtils::write10xCounts("./kb_out/counts_filtered", gene.symbol = genes[,2], filtered_m1, overwrite = T) # write result

cat("Filtering completed correctly.\n")
