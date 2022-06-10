FROM continuumio/miniconda3
# https://hub.docker.com/r/continuumio/miniconda3/tags
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
RUN rm -f sratoolkit.2.11.3-ubuntu64.tar.gz


# create conda environment
ADD pre_process.yml .
RUN conda init bash \
    && . ~/.bashrc \
    && conda env create -f pre_process.yml \
    && conda activate pre_process



ADD /ref/ /lib/ref/
ADD /scripts/* /usr/local/bin/
RUN chmod +x /usr/local/bin/*
CMD [ "/usr/local/bin/pre_process.sh", "/var/samples.csv"]
#RUN source /scripts/pre_process.sh