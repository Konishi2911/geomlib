FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y make wget gcc-10 g++-10 libgtest-dev libgmock-dev git

WORKDIR /tmp
RUN wget https://github.com/Kitware/CMake/releases/download/v3.21.1/cmake-3.21.1.tar.gz && tar xvfz cmake-3.21.1.tar.gz

WORKDIR /tmp/cmake-3.21.1
RUN export CC=gcc-10 && export CXX=g++-10 \
	&& ./bootstrap -- -DCMAKE_USE_OPENSSL=OFF && make && make install

RUN apt-get install -y libopenblas-dev liblapack-dev liblapacke-dev