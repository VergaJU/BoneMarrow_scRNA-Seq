FROM rocker/rstudio:4.1

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
    zlib1g-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    graphviz

RUN ls

RUN Rscript -e 'install.packages(c("rmarkdown", "tidyverse"), dependencies = TRUE, repos = "http://cran.us.r-project.org")'


RUN R --version
FROM bioconductor/bioconductor_docker:devel
RUN R --version

#RUN Rscript -e 'install.packages("Seurat"), dependencies = TRUE, repos = "http://cran.us.r-project.org")'
 
RUN Rscript -e 'BiocManager::install(c("slingshot", "SingleCellExperiment", "tradeSeq", "DelayedMatrixStats"))'