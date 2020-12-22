FROM ubuntu:20.04
RUN apt-get update \
 && apt-get install -y build-essential \ 
                       python3-dev \
                       python3-pip \
 && python3 -m pip install --upgrade setuptools 
                                     
COPY . /opt/UNCALLED/

RUN cd /opt/UNCALLED \
 && python3 setup.py install
