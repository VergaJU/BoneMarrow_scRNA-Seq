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
    graphviz

FROM bioconductor/bioconductor_docker:devel

RUN Rscript -e 'install.packages(c("rmarkdown", "markdown", "tidyverse", \
                "pheatmap", "Seurat",\
                "VennDiagram", "devtools"), dependencies = TRUE, repos = "http://cran.us.r-project.org")'


RUN Rscript -e 'BiocManager::install("AUCell")'
RUN Rscript -e 'BiocManager::install("RcisTarget")'
RUN Rscript -e 'BiocManager::install("GENIE3")'
RUN Rscript -e 'BiocManager::install("ComplexHeatmap")'

RUN Rscript -e 'install.packages("zoo")'
RUN Rscript -e 'install.packages("mixtools")'
RUN Rscript -e 'install.packages("rbokeh")'
RUN Rscript -e 'install.packages("DT")'
RUN Rscript -e 'install.packages("NMF")'
RUN Rscript -e 'install.packages("R2HTML")'

RUN Rscript -e 'install.packages("Rtsne")'
RUN Rscript -e 'install.packages("doMC")'
RUN Rscript -e 'install.packages("doRNG")'

RUN Rscript -e ' devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)' 
RUN Rscript -e ' devtools::install_github("aertslab/SCENIC")' 


WORKDIR /rstudio