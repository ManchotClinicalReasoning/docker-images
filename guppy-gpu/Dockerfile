FROM nvidia/cuda:10.1-runtime-ubuntu16.04

ARG PACKAGE_VERSION=3.4.4
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install --yes apt-transport-https \
                          libcurl4-openssl-dev \
                          libssl-dev \
                          libhdf5-cpp-11 \
                          libzmq5 \
                          libboost-atomic1.58.0 \
                          libboost-chrono1.58.0 \
                          libboost-date-time1.58.0 \
                          libboost-filesystem1.58.0 \
                          libboost-program-options1.58.0 \
                          libboost-regex1.58.0 \
                          libboost-system1.58.0 \
                          libboost-log1.58.0 \
                          libboost-iostreams1.58.0 \
                          wget

RUN    cd /tmp &&\
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_${PACKAGE_VERSION}-1~xenial_amd64.deb && \
    dpkg -i /tmp/ont_guppy_${PACKAGE_VERSION}-1~xenial_amd64.deb && \
    rm *.deb && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
