
FROM centos:7

RUN yum update -y

# Add packages for EH build
RUN yum install -y \
    bzip2-devel \
    git \
    libcurl-devel \
    libstdc++-static \
    xz-devel \
    zlib-devel

# Add packages for cmake build
RUN yum install -y \
    gcc \
    gcc-c++ \
    make \
    openssl \
    openssl-devel \
    wget

# Build newer cmake from source
RUN mkdir -p /download/cmake-3.20.5 && \
    cd /download/cmake-3.20.5 && \
    wget https://cmake.org/files/v3.20/cmake-3.20.5.tar.gz && \
    tar -xzf cmake-3.20.5.tar.gz && \
    mkdir build && cd build && \
    ../cmake-3.20.5/configure && \
    make -j8 && make install && \
    cd && rm -rf /download

# Build gcc9 from src
RUN yum install -y bzip2
RUN mkdir -p /download/gcc-9.4.0 && cd /download/gcc-9.4.0 && wget ftp://ftp.gnu.org/gnu/gcc/gcc-9.4.0/gcc-9.4.0.tar.xz && \
    tar -xJf gcc-9.4.0.tar.xz && cd gcc-9.4.0 && ./contrib/download_prerequisites && \
    cd .. && mkdir -p build && cd build && \
    ../gcc-9.4.0/configure \
        --prefix=/opt/gcc-9.4.0 \
        --disable-bootstrap \
        --disable-multilib \
        --enable-lto \
        --with-system-zlib \
        --enable-languages=c,c++ && \
    make -j8 && make install && cd / && rm -rf /download

# Set CC/CXX to new gcc version:
ENV CC /opt/gcc-9.4.0/bin/gcc
ENV CXX /opt/gcc-9.4.0/bin/g++

