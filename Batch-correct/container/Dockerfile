FROM nvidia/cuda:12.2.2-devel-ubuntu22.04

# Basic libraries

RUN apt-get update --fix-missing && apt-get install -y apt-utils\
    procps\
    wget\
    bzip2\
    ca-certificates\
    libglib2.0-0\
    libxext6\
    libsm6\
    libxrender1\
    git\
    mercurial\
    subversion\
    build-essential\
    libcurl4-gnutls-dev\
    libssl-dev\
    libxml2-dev\
    libhdf5-dev\
    graphviz

RUN apt -y install software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

RUN apt-get update && \
    apt-get install -y tzdata && \
    ln -fs /usr/share/zoneinfo/$TZ /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata


RUN apt -y install python3.9
RUN apt update && apt -y install curl
RUN curl -o get-pip.py "https://bootstrap.pypa.io/get-pip.py"
# RUN ln -sf /usr/bin/python3.9 /usr/bin/python3

RUN python3 get-pip.py
RUN pip install scanpy scvi-tools matplotlib seaborn torch
RUN pip install leidenalg 
RUN pip install scikit-misc
RUN pip install igraph
RUN pip install scib_metrics
RUN pip install optuna
RUN pip install plotly kaleido

RUN pip install scib
WORKDIR /usr/local/lib/python3.10/dist-packages/scib/knn_graph
RUN  g++ -std=c++11 -O3 knn_graph.cpp -o knn_graph.o
WORKDIR /