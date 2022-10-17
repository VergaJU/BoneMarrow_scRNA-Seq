FROM ubuntu:20.04
# Basic libraries
RUN apt-get update --fix-missing && apt-get install -y apt-utils procps wget bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion build-essential libcurl4-gnutls-dev libssl-dev libxml2-dev libhdf5-dev graphviz

# Miniconda
ENV PATH /opt/conda/bin:$PATH
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc


# Install python and required libraries
ENV PYTHON_VERSION=3.8.10
RUN conda install python=${PYTHON_VERSION}
ADD requirements.txt .
RUN pip3 install -r requirements.txt

# Install R packages

RUN Rscript -e 'install.packages(c("tidyverse", "scPOP", "optparse"))'

RUN chmod -R o+rX /opt/conda/lib/python3.8/site-packages/*
