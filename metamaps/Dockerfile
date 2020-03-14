FROM debian:stable-slim

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y autoconf build-essential git zlib1g-dev libboost-all-dev=1.67.0.1

#Create "virtual" boost library that matches the metaMaps build specification
RUN ln -s /usr/lib/x86_64-linux-gnu/libboost_math_c99.a /usr/lib/libboost_math_c99.a

RUN git clone https://github.com/DiltheyLab/MetaMaps.git && \
    cd MetaMaps && \
    git checkout 98102e9e684efa6a9903d8abe93600132c101ad0
    
WORKDIR /MetaMaps

RUN ./bootstrap.sh && \ 
    ./configure --with-boost=/usr && \ 
    make metamaps

#RUN apt-get -y autoremove build-essential autoconf

ENV PATH="/MetaMaps:${PATH}"

#Install Perl and dependencies for provided scripts

RUN apt-get install -y perl liblist-moreutils-perl libstatistics-basic-perl && \
    apt-get clean && \
    rm -rf /var/apt/lists
