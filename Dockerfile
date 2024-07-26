FROM ubuntu:20.04

# ensure local python is preferred over distribution python
ENV PATH /usr/local/bin:$PATH

# http://bugs.python.org/issue19846
# > At the moment, setting "LANG=C" on a Linux system *fundamentally breaks Python 3*, and that's not OK.
ENV LANG C.UTF-8

# Build-time environmental variable so that apt doesn't complain
ARG DEBIAN_FRONTEND=noninteractive

# dependencies
RUN apt update --fix-missing 
RUN apt install -y g++-10 libnlopt-cxx-dev bc python3.9 python3-pip python3-setuptools python3-distutils

# Install pybdei
RUN cd /usr/local/ && pip3 install --no-cache-dir numpy && pip3 install --no-cache-dir pybdei==0.8 && pip3 install --no-cache-dir pandas

# File Author / Maintainer
MAINTAINER Anna Zhukova <anna.zhukova@pasteur.fr>

# Clean up
RUN mkdir /pasteur

# The entrypoint runs BDEI parameter inference with command line arguments
ENTRYPOINT ["bdei_infer"]
