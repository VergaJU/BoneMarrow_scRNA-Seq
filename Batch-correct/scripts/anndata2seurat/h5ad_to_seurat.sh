#!/usr/bin/bash

file=$1 # input file, seurat object.rds


# Define paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"  # Get the script's directory
R_SCRIPT="$SCRIPT_DIR/seurat_to_mtx.R"
PY_SCRIPT="$SCRIPT_DIR/mtx_to_h5ad.py"

mkdir temp
docker run --group-add $(id -g) --rm -v $SCRIPT_DIR:/scripts -v $(pwd):/data:rw -w /data vergaju/bm_py_env python /scripts/anndata210x.py ${file}

docker run --group-add $(id -g) --rm -v $SCRIPT_DIR:/scripts -v $(pwd):/data -w /data vergaju/bm_r_env Rscript /scripts/10x2seurat.R ${file%h5ad}Rds # convert to 10 matrix
rm -fr temp