FROM debian:stable-slim
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y autoconf build-essential perl zlib1g-dev

COPY applyMetaMaps.pl /scripts/
ADD metamaps.tar.gz /
ADD boost_1_60_0.tar.gz /home/

# Install Boost
WORKDIR /home/boost_1_60_0
RUN ./bootstrap.sh --prefix=/usr/local --with-libraries=program_options
RUN ./b2 install
WORKDIR /home
RUN rm -rf boost_1_60_0

WORKDIR /metaMaps

RUN ./bootstrap.sh
RUN ./configure --with-boost=/usr/local/lib

WORKDIR /scripts
