#!/home/jacopo/anaconda3/envs/harmony/bin/Rscript
library(Seurat)
library(harmony)
library(future)
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

plan("multicore")

dat <- readRDS(opt$input_file) %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>%
    RunPCA(features = VariableFeatures(object = dat), verbose = F, seed.use = 1)

dat1 <- dat
dat1 <- FindNeighbors(dat1, dims = 1:20)
dat1 <- RunUMAP(dat1, dims = 1:20, seed.use = 1)
dat1 <- FindClusters(dat1, resolution = 0.2, random.seed = 1, verbose = T)

plot <- DimPlot(dat1, reduction="umap", label=T, pt.size=.1,raster=F)
svg("UMAP_cluster.svg")
print(plot)
dev.off()

plot <- DimPlot(dat1, reduction="umap",group.by="batch", label=T, pt.size=.1,raster=F)
svg("UMAP_batch.svg")
print(plot)
dev.off()

dat <- dat %>% 
    RunHarmony("batch", plot_convergence = TRUE)

dat <- FindNeighbors(dat, reduction = "harmony", dims = 1:20)
dat <- RunUMAP(dat, reduction = "harmony", dims = 1:20, seed.use = 1)
dat <- FindClusters(dat, resolution = 0.2, random.seed = 1, verbose = T)

plot <- DimPlot(dat, reduction="umap", label=T, pt.size=.1,raster=F)
svg("UMAP_cluster_harmony.svg")
print(plot)
dev.off()

plot <- DimPlot(dat, reduction="umap",group.by="batch", label=T, pt.size=.1,raster=F)
svg("UMAP_batch_harmony.svg")
print(plot)
dev.off()

saveRDS(dat1, file = opt$output_file)