FROM ubuntu:latest

RUN apt update -y && \
    apt install -y wget bzip2 libboost-tools-dev

# RUN rm -rf gpkg/boost && \
#     rm -rf gpkg/libs && \
#     wget -P /tmp https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.bz2 && \
#     tar xjpf /tmp/boost_1_88_0.tar.bz2 -C /tmp && \
#     bcp --boost=/tmp/boost_1_88_0/ --scan boost.hpp ./gpkg/

# apt install -y cmake-curses-gui cmake make gcc g++ libc-dev