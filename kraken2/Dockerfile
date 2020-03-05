FROM debian:stable-slim

LABEL website="https://github.com/DerrickWood/kraken2"
LABEL license="https://github.com/DerrickWood/kraken2/blob/master/LICENSE"

RUN apt-get update && apt-get -y install \
  git \
  zlib1g-dev \
  make \
  g++

RUN git clone https://github.com/DerrickWood/kraken2.git && \
  cd kraken2 && \
  ./install_kraken2.sh .

ENV PATH="/kraken2/kraken2:/kraken2/kraken2-build:/kraken2/kraken2-inspect:${PATH}"
