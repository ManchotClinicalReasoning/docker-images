FROM debian:stable-slim

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y autoconf build-essential perl=5.28.1-6 liblist-moreutils-perl zlib1g-dev libboost-all-dev=1.67.0.1 && \
    apt-get clean && \
    rm -rf /var/apt/lists

#Create "virtual" boost library that matches the metaMaps build specification
RUN ln -s /usr/lib/x86_64-linux-gnu/libboost_math_c99.a /usr/lib/libboost_math_c99.a

ADD metamaps.tar.gz /

WORKDIR /metaMaps

RUN ./bootstrap.sh && \ 
    ./configure --with-boost=/usr && \ 
    make metamaps

#RUN apt-get -y autoremove build-essential autoconf

ENV PATH="/metaMaps:${PATH}"

COPY applyMetaMaps.pl /scripts/
