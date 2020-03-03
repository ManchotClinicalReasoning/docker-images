FROM debian:stable-slim
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y perl liblist-moreutils-perl r-base
COPY nanoporeRunStats_oneSample.pl /scripts/
COPY plotRunStats_oneSample.R /scripts/

