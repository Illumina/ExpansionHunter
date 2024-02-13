From ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

ENV GIT_SSL_NO_VERIFY=1

RUN apt-get -qq update && apt-get install -yq \
    autoconf \
    automake \
    build-essential \
    cmake \
    git-all \
    libboost-all-dev \
    libfreetype6-dev \
    liblzma-dev \
    libpng-dev \
    m4 \
    python3-pip \
    python-pkgconfig \
    software-properties-common \
    zlib1g-dev \
    libbz2-dev \
    googletest \
    ca-certificates && update-ca-certificates

ADD . /opt/expansionhunter

RUN mkdir /opt/expansionhunter/build
WORKDIR /opt/expansionhunter/build
RUN cmake /opt/expansionhunter -DCMAKE_INSTALL_PREFIX=/opt/expansionhunter && make && make install

CMD ["/opt/expansionhunter/build/ExpansionHunter"]