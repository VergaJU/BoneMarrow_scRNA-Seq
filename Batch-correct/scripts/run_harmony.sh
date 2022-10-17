#!/usr/bin/bash

exec > >(tee -ia harmony.LOG)

theta=(1 2 3)
lambda=(0.5 1 1.5)
sigma=(0.05 0.1 0.15)


for t in ${theta[@]};do 
    for l in ${lambda[@]};do 
        for s in ${sigma[@]};do 
            ../scripts/run_harmony.R -i DATASET_seurat_PCA.Rds -t ${t} -l ${l} -s ${s}
        done
    done
done
