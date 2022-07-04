#!/home/jacopo/anaconda3/envs/singler/bin/Rscript
library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(stringr)
library(BiocParallel)


multicoreParam <- MulticoreParam(workers = 6)
input_file <- commandArgs(trailingOnly = TRUE)[1]
reference <- commandArgs(trailingOnly = TRUE)[2]

dat <- readRDS(input_file)
ref <- readRDS(reference)

dat1 <- as.SingleCellExperiment(dat)
#ref <- as.SingleCellExperiment(ref)

# run SingleR
pred <- SingleR(test=dat, ref=ref, labels=ref[['celltype']], de.method="wilcox", assay.type.test = 'logcounts', assay.type.ref = 'logcounts', BPPARAM=multicoreParam)

df <- data.frame(label = pred$labels)
rownames(df) <- rownames(pred)

# save results

dat <- AddMetaData(dat, df)

#plot <- DimPlot(dat, reduction="umap",group.by="label", label=T, pt.size=.1,raster=T)
#plot <- plot + NoLegend()
#png(str_replace(dat, ".Rds", "_UMAP_label.png"))
#print(plot)
#dev.off()


saveRDS(pred, file = str_replace(input_file, ".Rds", "_singleR_results.Rds"))
saveRDS(dat, file = str_replace(input_file, ".Rds", "_labelled.Rds"))
write.table(df, file =str_replace(input_file, ".Rds", '_cell_labels.csv'), sep = '\t', row.names = T, col.names = T, quote = F)