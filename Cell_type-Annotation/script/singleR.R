#!/home/jacopo/anaconda3/envs/singler/bin/Rscript


library(optparse)


option_list = list(
  make_option(
    c("-r", "--ref_loc"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to reference files.'
  ),
  make_option(
    c("-s", "--sample_loc"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to sample file, the one to be annotated.'
  ),
  make_option(
    c("-l", "--label"),
    action = "store",
    default = "label",
    type = 'character',
    help = 'label to assig cell types'
  ),
  make_option(
    c("-c","--celltype_ref"),
    action="store",
    default = "ct",
    type = "character",
    help="label with reference cell types"
  ),
  make_option(
    c("-t","--threads"),
    action="store",
    default = 0,
    type="integer",
    help='Number of threads, default:0 (all).'
  )
)


opt <- parse_args(OptionParser(option_list=option_list))


suppressMessages({library(Seurat)
  library(SingleCellExperiment)
  library(SingleR)
  library(stringr)
  library(BiocParallel)
  library(parallel)
  library(tidyverse)
})

if(opt$t==0){
  multicoreParam <- MulticoreParam(workers = detectCores())
} else{
  multicoreParam <- MulticoreParam(workers=opt$t)
}

# ref_loc <- commandArgs(trailingOnly = TRUE)[1]
# sample_loc <- commandArgs(trailingOnly = TRUE)[2]
# label <- commandArgs(trailingOnly = TRUE)[3]
sample_id <- str_extract(str_extract(opt$sample_loc, '[^/]+$'), '[^\\.]+')
print(paste('Working on', sample_id))

dat <- readRDS(opt$sample_loc)
ref <- readRDS(opt$ref_loc)

dat1 <- as.SingleCellExperiment(dat)
ref <- as.SingleCellExperiment(ref)

# run SingleR
pred <- SingleR(test=dat1, ref=ref, labels=ref[[opt$celltype_ref]], de.method="wilcox", assay.type.test = 'logcounts', assay.type.ref = 'logcounts', BPPARAM=multicoreParam)

# Props polot

data_label <- as.data.frame(table(pred$labels)) %>%
  mutate(prop = Freq/sum(Freq)) %>%
  arrange(-prop)


file_path=tools::file_path_sans_ext(opt$sample_loc)

ggplot(data = data_label, aes(x=Var1,y=prop)) +
  geom_bar(stat='identity') + 
  theme(axis.text.x = element_text(angle=45,vjust = 1, hjust = 1),
        plot.margin = margin(1,1,1,1, "cm"))+
  labs(title="Frequency cell types",
       x='cell type',
       y='Proportion') + 
  scale_x_discrete(limits = data_label$Var1)


ggsave(paste0(file_path,"_proportion_celltypes.png"),width = 10, height = 6)


df <- data.frame(label = pred$labels)
rownames(df) <- rownames(pred)
colnames(df) <- opt$label
# save results

dat <- AddMetaData(dat, df)

#plot <- DimPlot(dat, reduction="umap",group.by="label", label=T, pt.size=.1,raster=F)
#plot <- plot + NoLegend()
#svg(str_replace(input_file, ".Rds", "_UMAP_label.svg"))
#print(plot)
#dev.off()


saveRDS(pred, file = paste0(file_path, "_singleR_results.Rds"))
saveRDS(dat, file = paste0(file_path, "_labelled.Rds"))
write.table(df, file =paste0(file_path,  '_cell_labels.csv'), sep = '\t', row.names = T, col.names = T, quote = F)