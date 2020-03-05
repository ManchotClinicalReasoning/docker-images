FROM debian:stable-slim
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y perl liblist-moreutils-perl r-base
COPY quick_AMR_analysis.pl /scripts/

