FROM biopython/biopython:latest

COPY scripts /scripts
RUN pip3 install --upgrade pip
RUN apt-get remove python3-numpy -y
RUN pip3 install numpy --upgrade
RUN pip3 install pandas --upgrade
RUN pip3 install scipy --upgrade
RUN pip3 install seaborn --upgrade

