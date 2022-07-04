library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(stringr)

dat <- readRDS(commandArgs(trailingOnly = TRUE)[1])
ref <- readRDS(commandArgs(trailingOnly = TRUE)[1])

dat1 <- as.SingleCellExperiment(dat)
ref <- as.SingleCellExperiment(ref)

# run SingleR
pred <- SingleR(test=dat1, ref=ref, labels=ref[['celltype']], de.method="wilcox", assay.type.test = 'logcounts', assay.type.ref = 'logcounts')

df <- data.frame(label = pred$labels)
rownames(df) <- rownames(pred)

# save results

dat <- AddMetaData(dat, df)

plot <- DimPlot(dat, reduction="umap",group.by="label", label=T, pt.size=.1,raster=F)
plot <- plot + NoLegend()
svg(str_replace(input_file, ".Rds", "_UMAP_label.svg"))
print(plot)
dev.off()


saveRDS(pred, file = str_replace(input_file, ".Rds", "_singleR_results.Rds"))
saveRDS(dat, file = str_replace(input_file, ".Rds", "_labelled.Rds"))
write.table(df, file =str_replace(input_file, ".Rds", '_cell_labels.csv'), sep = '\t', row.names = T, col.names = T, quote = F)
