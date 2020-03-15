FROM debian:stable-slim

LABEL website="https://github.com/lh3/minimap2 and https://github.com/NCBI-Hackathons/NovoGraph"
LABEL license="https://github.com/lh3/minimap2/blob/master/LICENSE.txt and https://github.com/NCBI-Hackathons/NovoGraph/blob/master/LICENSE"


RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y perl liblist-moreutils-perl libstatistics-basic-perl r-base git


#Install minimap2
RUN git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && \
    git checkout tags/v2.17

RUN cd minimap2 && make

ENV PATH="/minimap2:${PATH}"


#Install (clone) novograph
RUN git clone https://github.com/NCBI-Hackathons/NovoGraph.git
RUN cd NovoGraph && \
    git checkout 3804a0f6a9d42e2b588052262a58b07f7412338d
RUN cd NovoGraph/src && \
    make all
ENV PATH="/NovoGraph/src:${PATH}"

COPY amr_1.pl /scripts/
COPY perlLib /scripts/perlLib

