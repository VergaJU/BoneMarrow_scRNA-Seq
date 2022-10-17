#!/usr/bin/bash

export PATH="$PATH:/home/jacopo/Documents/GitHub/BoneMarrow_scRNA-Seq/Merge_bath-correct/scripts/"

exec > >(tee -ia harmony.LOG)

batch=(1 2 3)
cell=(0 1 2)


for b in ${batch[@]};do 
    for c in ${cell[@]};do 
        run_harmony_celltype.R -i ${1} -t ${b} -l ${c}
    done
done