FROM debian:stable-slim

LABEL website="https://github.com/DerrickWood/kraken2 and https://github.com/jenniferlu717/Bracken"
LABEL license="https://github.com/DerrickWood/kraken2/blob/master/LICENSE and https://github.com/jenniferlu717/Bracken/blob/master/LICENSE"

RUN apt-get update && apt-get -y install \
  git \
  g++ \
  make \
  python \
  rsync \
  zlib1g-dev 

RUN git clone https://github.com/DerrickWood/kraken2.git && \
  cd kraken2 && \
  git checkout 066a40c8f43f3b7fe9ac31a7e169baaaa610e5b9 && \
  ./install_kraken2.sh .

ENV PATH="/kraken2:${PATH}"

RUN git clone https://github.com/jenniferlu717/Bracken.git && \
  cd Bracken && \
  git checkout c3d499d6ef649c44b08c098f2babeb11c34c40fa && \
  sh install_bracken.sh

ENV PATH="/Bracken:${PATH}"

