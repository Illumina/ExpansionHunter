FROM ubuntu:16.04

RUN apt-get -qq update && apt-get install -yq \
    wget curl software-properties-common && \
    wget -O llvm.key https://apt.llvm.org/llvm-snapshot.gpg.key && apt-key add llvm.key && rm -f llvm.key

RUN apt-add-repository "deb http://apt.llvm.org/xenial/ llvm-toolchain-xenial-5.0 main"

RUN apt-get -qq update && apt-get install -yq \
  build-essential \
  cmake \
  zlib1g-dev \
  libbz2-dev \
  valgrind \
  cppcheck \
  clang-5.0 \
  clang-format-5.0 \
  clang-tidy-5.0 \
  libboost-all-dev \
  &&  \
  apt-get clean -y && \
  rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-5.0 1000 && \
    update-alternatives --install /usr/bin/clang clang /usr/bin/clang-5.0 1000 && \
    update-alternatives --config clang && \
    update-alternatives --config clang++

RUN apt-get -qq update && apt-get install -yq \
  git \
  &&  \
  apt-get clean -y && \
  rm -rf /var/lib/apt/lists/*
