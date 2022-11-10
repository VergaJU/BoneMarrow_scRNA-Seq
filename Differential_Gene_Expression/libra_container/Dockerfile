FROM rocker/rstudio:4.0.5

RUN apt-get update --fix-missing && apt-get install -y apt-utils \
    procps \
    wget \
    bzip2 \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    subversion \
    build-essential \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    graphviz


RUN Rscript -e 'install.packages(c("rmarkdown", "tidyverse", "BiocManager", "devtools", "Seurat"), dependencies = TRUE, repos = "http://cran.us.r-project.org")'

RUN Rscript -e 'BiocManager::install(c("edgeR", "DESeq2", "limma"))'

Run Rscript -e 'devtools::install_github("neurorestore/Libra")'