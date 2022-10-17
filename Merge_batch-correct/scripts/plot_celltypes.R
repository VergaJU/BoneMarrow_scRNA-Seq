# Load the libraries (from Sarah script + biomart)
library(tidyverse) # packages for data wrangling, visualization etc
library(Seurat) # scRNA-Seq analysis package
library(clustree) # plot of clustering tree 
library(ggsignif) # Enrich your 'ggplots' with group-wise comparisons
library(clusterProfiler) #The package implements methods to analyze and visualize functional profiles of gene and gene clusters.
library(org.Hs.eg.db) # Human annotation package neede for clusterProfiler
library(ggrepel) # extra geoms for ggplo2
library(patchwork) #multiplots
library(reticulate)


dat <- readRDS("COMPLETE_harmony_labelled.Rds")

dat

colors = c()
for (f in 1:24){
  colors[f] = "grey"
}

for (lab in 1:24){
  newcolors = colors
  newcolors[lab] = "blue"
  plot <- DimPlot(dat, reduction="umap", group.by = "label", cols = newcolors, label=T, pt.size=.1, raster=F)
  ggsave(plot, width = 12, height = 6, dpi = 300, filename = paste0(lab, "_plot.png"))
  newcolors = colors
}


plot <- DimPlot(dat, reduction="umap", group.by = "label", split.by = FetchData(dat, slot = "label"), label=T, pt.size=.1,raster=F)
plot <- plot + NoLegend()
plot

svg(str_replace(input_file, ".Rds", "_UMAP_label.svg"))
print(plot)
dev.off()
