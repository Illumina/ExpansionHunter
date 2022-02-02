
FROM fedora:34 

RUN yum update -y

# Add packages for EH build
RUN yum install -y \
    bzip2-devel \
    cmake \
    gcc \
    gcc-c++ \
    git \
    libcurl-devel \
    libstdc++-static \
    openssl-devel \
    xz-devel \
    zlib-devel

# Add extra debug packages
RUN yum install -y time 

