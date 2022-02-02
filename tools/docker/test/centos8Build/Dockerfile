
FROM centos:8

# Add std packages for EH build
RUN yum update -y && yum install -y \
    bzip2-devel \
    cmake \
    gcc \
    gcc-c++ \
    libcurl-devel \
    make \
    openssl-devel \
    xz-devel \
    zlib-devel

# Configure yum to find libstdc++-static and install
RUN yum update -y && \
    yum install -y dnf-plugins-core && \
    yum config-manager --set-enabled powertools && \
    yum install -y libstdc++-static

