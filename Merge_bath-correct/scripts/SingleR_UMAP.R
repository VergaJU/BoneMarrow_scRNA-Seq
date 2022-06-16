library(Seurat)
library(SingleCellExperiment)
library(SingleR)

dat <- readRDS(opt$input_file)
ref <- readRDS(opt$ref)

dat <- as.SingleCellExperiment(dat)
ref <- as.SingleCellExperiment(ref)

# run SingleR
pred <- SingleR(test=dat, ref=ref, labels=ref[['celltype']], de.method="wilcox", assay.type.test = 'logcounts', assay.type.ref = 'logcounts')

df <- data.frame(label = pred$labels)
rownames(df) <- rownames(pred)

# save results
pred <- as.Seurat(pred)


plot <- DimPlot(pred, reduction="umap",group.by="celltype", label=T, pt.size=.1,raster=F)
plot <- plot + NoLegend()
svg("UMAP_batch_harmony.svg")
print(plot)
dev.off()


saveRDS(pred, file = paste0(opt$output_name, ".Rds"))
write.table(df, file = paste0(opt$output_name, '_cell_labels.csv'), sep = '\t', row.names = T, col.names = T, quote = F)
