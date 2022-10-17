# Merge data

Merge all the files inputted from a text file containing the paths to be added. The file names are added in the variable `-b`. Output are a Seurat and a SingleCellExperiment files.

## Usage:

```
Usage: ./scripts/merge_data.R [options]


Options:
        -i INPUT_FILE, --input_file=INPUT_FILE
                Path to input file, it must be a file containing the paths for the files to be merged.

        -o OUTPUT_PREFIX, --output_prefix=OUTPUT_PREFIX
                Prefix for naming output file.

        -b BATCH_KEY, --batch_key=BATCH_KEY
                Variable where to store the batch names.

        -h, --help
                Show this help message and exit

```