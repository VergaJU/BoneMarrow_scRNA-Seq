#!/usr/bin/bash

file=$1 # input file, seurat object.rds

mkdir temp # create temporary directory

/home/jacopo/Documents/GitHub/BoneMarrow_scRNA-Seq/Batch-correct/scripts/seurat2anndata/seurat_to_mtx.R $file  # convert to 10x matrix

mv ./temp/genes.tsv ./temp/old.tsv # mv genes in a new file

# add the genes in a new genes.tsv file with 2 columns

for f in $(cat ./temp/old.tsv);do
    echo -e "$f\t$f" >> ./temp/genes.tsv
    done

#convert in a h5ad file
/home/jacopo/Documents/GitHub/BoneMarrow_scRNA-Seq/Batch-correct/scripts/seurat2anndata/mtx_to_h5ad.py $file

# remove temp directory
rm -fr temp/