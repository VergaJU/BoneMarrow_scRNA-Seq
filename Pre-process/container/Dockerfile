FROM r-base:4.0.0
# Basic libraries
RUN apt-get update --fix-missing && apt-get install -y apt-utils procps wget bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion build-essential libcurl4-gnutls-dev libssl-dev libxml2-dev libhdf5-dev graphviz


# Download and add in the path CellRanger
RUN wget -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.1.2.tar.gz?Expires=1654317776&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NTQzMTc3NzZ9fX1dfQ__&Signature=fpK4V4HtqhjFljx4437ou0E9Lzsfmzhfa88-p-GQZiJDWUU4H5fKyIEPAXmvRnBY1e09IKjk-CqwLxV~Xi3c1Q71Bx8EIIFfz7FW~s0Yjxp-as-UkoYOsOtHAXdSbosWCVM1JtWLMnUVADIYuN5nsCct4wkpyEmiPm1woQU~vlQL28Mej-TlcHfpFgpxZHpDJ5CUQMQ804tJ1HivwF1mEoaF4Vpf3MCQAYvcZ-JAy32rz3zV06u7Nd3i5N6c9YV3HRZ7u5CgoyrYe7NwAFgJ4dXfgG-5VUAd6eKgIgew5UsW8celvwc8HHCln8m3Um5cRrkHMixoHH71z1babKo7KQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
RUN tar -zxvf cellranger-6.1.2.tar.gz

ENV PATH=$HOME/cellranger-6.1.2:$PATH

#Download and install SraToolkit
RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.4.1/sratoolkit.2.4.1-ubuntu64.tar.gz
# Unzip the archive
RUN tar xzvf sratoolkit.2.4.1-ubuntu64.tar.gz
ENV PATH=$HOME/directory/sratoolkit.2.4.1-ubuntu64/bin:$PATH

# Miniconda
ENV PATH /opt/conda/bin:$PATH
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ADD pre_process.yml .
RUN conda env create -f pre-process.yml