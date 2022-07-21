# Pre processing

The pre processing consist on:
- Data retrieval:
    - sratoolkit fasterq-dump for sra entries
    - wget for bam entries
- Conversion of the files in fastq format:
    - sratoolkit fasterq-dump for sra entries
    - cellranger for bam files
- Pseudoalignment and counts using kallisto|bustools
- Empty droplets filtration with EmptyDrops
- Doublets removal using scDblFinder

## Usage: 

All the pipeline is organised in a Docker container to ensure reproducibility and portability of it, the container is available [here](https://hub.docker.com/repository/docker/vergaju/pre-process).

To obtain the latest version of it run:
```
docker push vergaju/pre-process:v5
```

The container contains the seguent softwares:
- kallisto|bustools python wrapper (0.27.0)
- sratoolkit (2.11.3)
- cellranger (6.1.2)
- DropletUtils (1.14.0)
- scDblFinder (1.8.0)

The index is the built-in from kallisto (from Jan 2022).

The usage is:
```
docker run pre-process <csv file>
```

The input file is a csv files organised as follow:
- patient or experiment
- name of the run
- entry of the run/link
- technology

example:

```
SAMN18822752,SRR14295357,SRR14295357,10xv3
SAMN18822743,SRR14295358,<link to bam file>,10xv3
```

Headers aren't needed and the technology supported at today are:
- 10xv2
- 10xv3

## Steps:

For each row (entry) two nested folders will be created:
```
Parent Dir
└── Patient
     └── Name
```

### Data retrival

Bam files created with CellRanger have to be inserted as links, wget will download the file. CellRanger then will convert back the bam file to fastq files:

```
# download
wget -c <link>
filename=$(ls) # update filename

## Use cellranger to get the fastq files
~/.cellranger-6.1.2/lib/bin/bamtofastq \
    --nthreads=6 \
    --reads-per-fastq=100000000000 \
    --traceback \
    ${filename} \
    ./fastq/ \
    |& tee -a bamtofastq.LOG
```

The log file will be saved in the file `bamtofastq.LOG`.

Sra entries will be automatically downloaded and converted using `fasterq-dump`:

```
fasterq-dump <entry> -o <entry> -O fastq/ -S --include-technical |& tee -a fasterq-dump.LOG
cache-mgr -c ~/.sratoolkit.2.11.3-ubuntu64/cache/
```

The fastq files will be named as the entry. The log file frm fasterq-dump will be saved as `fasterq-dump.LOG`.

### Pesudoalignment and counts:

All the fastq file will be processed with kallisto|bustools:

```
kb count --verbose -t 6 \
    --cellranger \
    -i $index \
    -g $t2g \
    -x $tech \
    -o ./kb_out ${fastqs[@]} \
    --overwrite |& tee -a kallisto.LOG 
```

As previously the log file will be saved in `kallisto.LOG`.

### Empty/doublets filtering

Finally, the matrix obtained from kb fill be filtered from empty cells (FDR 0.1)  and doublets:

```
filter_empty_v2.R ./kb_out/counts_unfiltered/cellranger/ 0.1 |& tee -a filter_empty.LOG
```

the log file will be saved here: `filter_empty.LOG`

### Cleanup files

Since fastq and bus files are big and occupy a lot of storage, the last step is to remove them:

```
find . -name "*.fastq" -delete
find . -name "*.bus" -delete
```

The final structure of the directory for each row will be as follow:
```
.
├── main.LOG
└── SAMN15892704
    └── SRR12506863
        ├── fasterq-dump.LOG
        ├── fastq
        ├── filter_empty.LOG
        ├── kallisto.LOG
        └── kb_out
              ├── 10x_version3_whitelist.txt
              ├── counts_filtered
              │   ├── barcodes.tsv
              │   ├── genes.tsv
              │   └── matrix.mtx
              ├── counts_unfiltered
              │   ├── cellranger
              │   │   ├── barcodes.tsv
              │   │   ├── genes.tsv
              │   │   └── matrix.mtx
              │   ├── cells_x_genes.barcodes.txt
              │   ├── cells_x_genes.genes.txt
              │   └── cells_x_genes.mtx
              ├── inspect.json
              ├── kb_info.json
              ├── matrix.ec
              ├── run_info.json
              └── transcripts.txt
```

- `main.LOG` Contains the LOG file from the main script ([pre_process.sh](./scripts/pre_process.sh)), one for all the entries inputted.
- Each entry has it's log files inside the entry folder.
- `kb_out` is the folder containing all the results from the pre-processing:
    - `counts_unfiltered` contains the matrices from `kb count` in kallisto and CellRanger outputs
    - `counts_filtered` contains matrix, barcodes and genes filtered with EmptyDroplets and scDblFinder. These files will be considered for the next Quality Control step.


```mermaid
    flowchart TD;
      A[Create Entry folder]-->B[Read entries (BAM or SRA)];
      B-->C[is BAM];
      C-->D[Download with wget];
      B-->E[is SAM];
      E-->F[Download and convert with fasterq-dump];

```