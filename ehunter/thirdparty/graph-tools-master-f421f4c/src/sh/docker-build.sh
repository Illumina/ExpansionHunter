#!/usr/bin/env bash

# Make a build using Docker
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

set -e

PULL=1
CLANG=0
CLANG_ASAN=0
CLANG_MSAN=0
VALGRIND=0
while test $# -gt 0
do
    case "$1" in
        --no-pull) echo "Not pulling image."; PULL=0
            ;;
        --valgrind) echo "Using valgrind"; VALGRIND=1
            ;;
        --clang) echo "Using CLang"; CLANG=1
            ;;
        --clang-asan) echo "Using CLang+Address Sanitizer"; CLANG_ASAN=1; CLANG=1
            ;;
        --clang-msan) echo "Using CLang+Memory Sanitizer"; CLANG_MSAN=1; CLANG=1
            ;;
        *) echo "Unknown argument $1"; exit 1
            ;;
    esac
    shift
done

if [[ ${PULL} -ne 0 ]]; then
    docker pull ilmncgrpmi/cpp-dev:latest
fi

EXTRA_CMAKE=""
if [[ ${CLANG} -ne 0 ]]; then
    EXTRA_CMAKE="-DCMAKE_CXX_COMPILER=clang++"
fi
if [[ ${CLANG_ASAN} -ne 0 ]]; then
    EXTRA_CMAKE="${EXTRA_CMAKE} -DUSE_ASAN=ON"
fi
if [[ ${CLANG_MSAN} -ne 0 ]]; then
    EXTRA_CMAKE="${EXTRA_CMAKE} -DUSE_MSAN=ON"
fi

if [[ ${VALGRIND} -ne 0 ]]; then
    VALGRIND='for x in `find /graphtools-build/tests -executable -type f`; do valgrind --leak-check=full --xml=yes --xml-file=/graphtools-source/valgrind_`basename ${x}`.xml $x ; done && find /graphtools-source/valgrind_*.xml | xargs -n1 python /graphtools-source/src/sh/valgrind-check.py &&'
else
    VALGRIND=""
fi

docker run -v ${DIR}/../..:/graphtools-source --rm ilmncgrpmi/cpp-dev:latest /bin/bash -c "mkdir -p /graphtools-build &&  \
   cd /graphtools-build &&  \
   cmake ../graphtools-source -DBUILD_TESTS=ON -DCMAKE_INSTALL_PREFIX=/graphtools-install ${EXTRA_CMAKE} &&  \
   make -j8 && \
   make test && \
   make install && \
   ${VALGRIND} \
   cd /graphtools-install && \
   tar czf /graphtools-source/graphtools-install.tar.gz * "

