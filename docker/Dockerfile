FROM ubuntu:20.04

ARG git_checkout=master
ARG git_url=https://github.com/skovaka/UNCALLED.git
ARG uncalled_dir=/opt/UNCALLED

RUN apt-get update \
 && apt-get install -y build-essential \
                       git \ 
                       python3-dev \
                       python3-pip \
 && python3 -m pip install --upgrade setuptools 
                                     
RUN git clone --recursive $git_url $uncalled_dir \
 && cd $uncalled_dir \
 && git checkout $git_checkout \
 && python3 setup.py install
