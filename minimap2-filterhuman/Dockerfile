FROM debian:stable-slim

LABEL website="https://github.com/lh3/minimap2"
LABEL license="https://github.com/lh3/minimap2/blob/master/LICENSE.txt"


RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y perl liblist-moreutils-perl libstatistics-basic-perl r-base git


#Install minimap2
RUN git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && \
    git checkout tags/v2.17

RUN cd minimap2 && make

ENV PATH="/minimap2:${PATH}"

COPY removeHumanReads.pl /scripts/
COPY perlLib /scripts/perlLib

