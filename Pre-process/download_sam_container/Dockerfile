FROM ubuntu:20.04
# https://hub.docker.com/r/continuumio/miniconda3/tags
# Basic libraries
RUN apt-get update --fix-missing && apt-get install -y apt-utils procps wget bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion build-essential libcurl4-gnutls-dev libssl-dev libxml2-dev libhdf5-dev graphviz
#Download and install SraToolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
# Unzip the archive
RUN tar xzvf sratoolkit.2.11.3-ubuntu64.tar.gz
ENV PATH=/sratoolkit.2.11.3-ubuntu64/bin:$PATH
RUN vdb-config --restore-defaults
RUN rm -f /root/.ncbi/user-settings.mkfg
ADD user-settings.mkfg /root/.ncbi/user-settings.mkfg
RUN rm -f sratoolkit.2.11.3-ubuntu64.tar.gz

## Install samtools:

RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
RUN tar -xvf samtools-1.16.1.tar.bz2
RUN cd  samtools-1.16.1
RUN ./configure --prefix=/opt/
RUN make
RUN make install
RUN samtools -v

