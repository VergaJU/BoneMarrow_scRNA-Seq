#!/usr/bin/bash

file=$1 # input file, seurat object.rds


# Define paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"  # Get the script's directory
R_SCRIPT="$SCRIPT_DIR/seurat_to_mtx.R"
PY_SCRIPT="$SCRIPT_DIR/mtx_to_h5ad.py"

mkdir temp # create temporary directory

docker run --rm -v $SCRIPT_DIR:/scripts -v $(pwd):/data -w /data vergaju/bm_r_env Rscript /scripts/seurat_to_mtx.R ${file} # convert to 10 matrix

docker run --rm -v $SCRIPT_DIR:/scripts -v $(pwd):/data -w /data vergaju/bm_py_env python /scripts/mtx_to_h5ad.py ${file}

rm -fr temp