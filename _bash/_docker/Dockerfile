############################################################
# Dockerfile for Tanytarsus gracilentus genomics.
############################################################


FROM condaforge/mambaforge:24.9.2-0
MAINTAINER Lucas Nell <lucnell@gmail.com>
WORKDIR /app
ENV HOME /app
ARG DEBIAN_FRONTEND=noninteractive



SHELL ["/bin/bash", "--login", "-c"]

# --------------------------------
# Making everything accessible on the cluster where they make you non-root
# --------------------------------
# I want to always be able to access these:
RUN chmod -R 777 /opt/conda && \
    chmod -R 777 /app


# Basic requirements
RUN apt-get clean && \
    apt-get update -y && \
    apt-get install -y build-essential && \
    apt-get install -y locales wget rsync unzip nano sudo ocaml-nox zlib1g-dev \
        libbz2-dev liblzma-dev gsl-bin libgsl0-dev libboost-all-dev && \
    rm -rf /var/lib/apt/lists/*

# Set locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8
ENV TZ America/New_York

# simple python script for sliding window of coverage based on mpileup and
# some helper functions in bash:
COPY window-mpileup.py .
COPY helpers.sh .
RUN chmod +x window-mpileup.py && \
    mv window-mpileup.py /usr/local/bin/ && \
    echo -e "\n\n" >> ~/.bashrc && \
    cat helpers.sh >> ~/.bashrc && \
    echo -e "\n\n" >> ~/.bashrc && \
    rm helpers.sh

# Make primary and quality-control conda environments
COPY main-env.yml .
RUN mamba env create -q -f main-env.yml && \
    rm main-env.yml && \
    conda clean --all -y
COPY qc-env.yml .
RUN mamba env create -q -f qc-env.yml && \
    rm qc-env.yml && \
    conda clean --all -y


RUN conda init bash

# --------------------------------
# Manually install SNAPE-pooled
# --------------------------------
RUN cd /opt && \
	sha="284bfe054da37232d3e23dfbfa6f67bcaa9cf89a" && \
    wget -q https://github.com/EmanueleRaineri/snape-pooled/archive/${sha}.zip && \
	unzip -q ${sha}.zip && \
    rm ${sha}.zip && \
	mv snape-pooled-${sha} snape-pooled && \
	cd snape-pooled && \
	make && \
	mv snape-pooled /usr/local/bin/ && \
	cd .. && \
	rm -r snape-pooled


# --------------------------------
# Manually install npstat
# --------------------------------
RUN cd /opt && \
	sha="dc230fba3ac434ed6526ab7cd63a67d69d34b7d6" && \
    wget -q https://github.com/lucaferretti/npstat/archive/${sha}.zip && \
	unzip -q ${sha}.zip && \
    rm ${sha}.zip && \
	mv npstat-${sha} npstat && \
	cd npstat && \
	make && \
	mv npstat /usr/local/bin/ && \
	cd .. && \
	rm -r npstat
