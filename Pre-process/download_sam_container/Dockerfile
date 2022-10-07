FROM ubuntu:20.04

# don't freeze pls
ENV TZ=Europe/Dublin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Basic libraries
RUN apt-get update --fix-missing && apt-get install -y\
    apt-utils\
    wget\
    build-essential\
    libncurses5-dev\
    zlib1g-dev\
    librust-bzip2-dev\
    liblzma-dev\
    graphviz

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
WORKDIR /samtools-1.16.1/
RUN ./configure --prefix=/opt/
RUN make
RUN make install
WORKDIR /
ENV PATH=/opt/bin:$PATH