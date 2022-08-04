#!/usr/bin/bash

export PATH="$PATH:/home/jacopo/Documents/GitHub/BoneMarrow_scRNA-Seq/Merge_bath-correct/scripts/"

batch=(1 2 3)
cell=(0 1 2)


for b in ${batch[@]};do 
    for c in ${cell[@]};do 
        run_harmony_celltype.R -i DATASET_seurat_PCA.Rds -t ${b} -l ${c} |& tee -a harmony.LOG
    done
done
