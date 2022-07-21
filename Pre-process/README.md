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
- entry of the run
- technology

example:

```
SAMN18822752,SRR14295357,SRR14295357,10xv3
SAMN18822743,SRR14295358,SRR14295358,10xv3
SAMN18822742,SRR14295359,SRR14295359,10xv3
SAMN18822751,SRR14295360,SRR14295360,10xv3
SAMN18822750,SRR14295361,SRR14295361,10xv3
SAMN18822749,SRR14295362,SRR14295362,10xv3
SAMN18822748,SRR14295363,SRR14295363,10xv3
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