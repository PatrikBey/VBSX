FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

ENV LANG="C.UTF-8" \
    LC_ALL="C.UTF-8"

# Set up the environment
ENV OS=Linux 



# # update python3 for use in QC
RUN apt-get update
RUN apt-get -y install python3-pip
RUN pip3 install numpy==1.21 && \
    pip3 install werkzeug==1.0.1 && \
    pip3 install progressbar && \
    pip3 install scipy && \
    pip3 install matplotlib && \
    pip3 install pyreadr

RUN apt-get update && apt-get install -y r-base

