FROM ubuntu:18.04 
# https://hub.docker.com/_/ubuntu
# Basic libraries
RUN apt-get update --fix-missing && apt-get install -y apt-utils procps wget bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion build-essential libcurl4-gnutls-dev libssl-dev libxml2-dev libhdf5-dev graphviz


# Download and add in the path CellRanger
ADD .cellranger-6.1.2 ./.cellranger-6.1.2/
ENV PATH=$HOME/.cellranger-6.1.2:$PATH

#Download and install SraToolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
# Unzip the archive
RUN tar xzvf sratoolkit.2.11.3-ubuntu64.tar.gz
ENV PATH=$HOME/sratoolkit.2.11.3-ubuntu64/bin:$PATH
RUN vdb-config --restore-defaults
RUN rm -f /root/.ncbi/user-settings.mkfg
ADD user-settings.mkfg /root/.ncbi/user-settings.mkfg

# Miniconda
ENV PATH /opt/conda/bin:$PATH
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ADD pre_process.yml .
RUN conda env create -f pre_process.yml