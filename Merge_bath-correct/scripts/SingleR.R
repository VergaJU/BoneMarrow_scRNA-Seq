#!/home/jacopo/miniconda3/envs/singler/bin/Rscript

library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(stringr)
library(BiocParallel)

multicoreParam <- MulticoreParam(workers = 6)

ref_loc <- commandArgs(trailingOnly = TRUE)[1]
sample_loc <- commandArgs(trailingOnly = TRUE)[2]
label <- commandArgs(trailingOnly = TRUE)[3]
sample_id <- str_extract(str_extract(sample_loc, '[^/]+$'), '[^\\.]+')
print(paste('Working on', sample_id))

dat <- readRDS(sample_loc)
ref <- readRDS(ref_loc)

dat1 <- as.SingleCellExperiment(dat)
ref <- as.SingleCellExperiment(ref)

# run SingleR
pred <- SingleR(test=dat1, ref=ref, labels=ref[[label]], de.method="wilcox", assay.type.test = 'logcounts', assay.type.ref = 'logcounts', BPPARAM=multicoreParam)

df <- data.frame(label = pred$labels)
rownames(df) <- rownames(pred)

# save results

dat <- AddMetaData(dat, df)

#plot <- DimPlot(dat, reduction="umap",group.by="label", label=T, pt.size=.1,raster=F)
#plot <- plot + NoLegend()
#svg(str_replace(input_file, ".Rds", "_UMAP_label.svg"))
#print(plot)
#dev.off()


saveRDS(pred, file = str_replace(sample_loc, ".Rds", "_singleR_results.Rds"))
saveRDS(dat, file = str_replace(sample_loc, ".Rds", "_labelled.Rds"))
write.table(df, file =str_replace(sample_loc, ".Rds", '_cell_labels.csv'), sep = '\t', row.names = T, col.names = T, quote = F)