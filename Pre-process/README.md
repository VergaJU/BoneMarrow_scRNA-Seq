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

All the pipeline is organised in a Docker container to ensure reproducibility and portability of it, the container is available [here](https://hub.docker.com/repository/docker/vergaju/pre-process).

To obtain the latest version of it run:
```
docker push vergaju/pre-process:v5
```

The container contains the seguent softwares:
- sratoolkit (2.11.3)
- cellranger (6.1.2)
- DropletUtils (1.14.0)
- scDblFinder (1.8.0)


The usage is:
```
docker run <input file>
```

The input file is a csv files organised as follow: