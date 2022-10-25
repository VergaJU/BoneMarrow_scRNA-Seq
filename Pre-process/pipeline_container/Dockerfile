FROM continuumio/miniconda3
# https://hub.docker.com/r/continuumio/miniconda3/tags
# Basic libraries
RUN apt-get update --fix-missing && apt-get install -y apt-utils \
    procps \
    wget \
    bzip2 \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    git \
    mercurial \
    subversion \
    build-essential \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    graphviz


# Download, install and add in the path CellRanger

RUN wget -O cellranger-7.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.1.tar.gz?Expires=1666668728&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjAuMS50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NjY2Njg3Mjh9fX1dfQ__&Signature=IhcRaL3qL6uz2zqvSHrNzavhyJkfTxZEbu-LJZnHgkkJV~Ha-fh3~TFbtKPn-oRGreacpQdSZW~8aZdAjBVTRh4D97orvJYdD~ASKMiNAL7982u-O1ppjRY69U8OWt3F8SAqLNIUw9nMLvlIjXUgcED1auX5YMv6qH58w1-fq9B0MRbs1Hw0tN~Pa9VueXlrqYRnGdoBwZkvvYF2EUI1NBgKDKndfREq4tRrLlm6vqYqq54WZNVBc~~H2v0aWFlHJQwxRpDeTT-Z5x1Yj1wzBXVmXJf6~KsyOng0P1UgHXcV3FB8S7HRjc8p5g5cA~haqSvvAmNxbZf4ogdAjI~Sjw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
RUN tar -xzvf cellranger-7.0.1.tar.gz
RUN rm -f cellranger-7.0.1.tar.gz
ENV PATH=/cellranger-7.0.1:$PATH
ENV PATH=/cellranger-7.0.1/lib/bin:$PATH


# Download reference transcriptome
#RUN wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
#RUN tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
#RUN rm -f refdata-gex-GRCh38-2020-A.tar.gz
#RUN ls -R refdata-gex-GRCh38-2020-A

#Download and install SraToolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
# Unzip the archive
RUN tar xzvf sratoolkit.2.11.3-ubuntu64.tar.gz
ENV PATH=/sratoolkit.2.11.3-ubuntu64/bin:$PATH
RUN vdb-config --restore-defaults
RUN rm -f /root/.ncbi/user-settings.mkfg
ADD user-settings.mkfg /root/.ncbi/user-settings.mkfg
RUN rm -f sratoolkit.2.11.3-ubuntu64.tar.gz


# create conda environment
RUN conda update --all
ADD pre_process.yml .
RUN conda init bash \
    && . ~/.bashrc \
    && conda env create -f pre_process.yml \
    && conda activate pre_process



#ADD /scripts/* /usr/local/bin/
#RUN chmod +x /usr/local/bin/*
#CMD [ "/usr/local/bin/pre_process.sh", "/var/samples.csv"] don't run if using with singularity
#RUN source /scripts/pre_process.sh