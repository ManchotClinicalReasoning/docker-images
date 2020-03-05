FROM debian:stable-slim

LABEL website="https://github.com/fenderglass/Flye"
LABEL license="https://github.com/fenderglass/Flye/blob/flye/LICENSE"

RUN apt-get update && apt-get -y install \
  gcc \
  git \
  g++ \
  make \
  python3 \
  python3-setuptools \
  zlib1g-dev


RUN git clone https://github.com/fenderglass/Flye.git && \
  cd Flye && \
  python3 setup.py install
