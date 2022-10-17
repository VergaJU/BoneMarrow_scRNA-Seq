#!/usr/bin/bash

file=$1 # input file, seurat object.rds

mkdir temp # create temporary directory

~/Documents/PROJECT/scripts_env/scripts/seurat_to_mtx.R $file  # convert to 10x matrix

mv ./temp/genes.tsv ./temp/old.tsv # mv genes in a new file

# add the genes in a new genes.tsv file with 2 columns

for f in $(cat ./temp/old.tsv);do
    echo -e "$f\t$f" >> ./temp/genes.tsv
    done

#convert in a h5ad file
~/Documents/PROJECT/scripts_env/scripts/mtx_to_h5ad.py $file

# remove temp directory
rm -fr temp/