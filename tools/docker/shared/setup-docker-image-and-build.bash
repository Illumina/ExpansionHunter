#!/usr/bin/env bash

set -o errexit
set -o nounset


if [ "$#" -ne 2 ]; then
    echo "ERROR $0 requires 2 arguments for minDockerGiB baseScript"
    exit 1
fi

minDockerGiB=$1
baseScriptDir="$(cd "$(dirname $2)" && pwd -P)"

tag1="eh$(basename "$(dirname "$baseScriptDir")")"
tag2="$(basename "$baseScriptDir")"

if [ "$EUID" -ne 0 ]; then
    echo "This script must be run with sudo"
    exit 2
fi

dockermem=$(docker info | awk '/Total Memory/ {print $3}')

if [ -z "$dockermem" ]; then
    echo "ERROR: No docker daemon found"
    exit 1
fi

if echo $dockermem | awk -v mingb="${minDockerGiB}" -F "GiB" '$1~/^[0-9.]+$/ && $1>=mingb {exit 1}' ; then
    echo "ERROR: Docker is configured with $dockermem of total memory. This build may be unstable when less than ${minDockerGiB}GiB is available."
    exit 1
fi

scriptDir=$(cd "$(dirname $0)" && pwd -P)
rootDir=$(cd $scriptDir && cd ../../.. && pwd -P)
dockerRoot=/builder
tag=${tag1}:${tag2}
buildDirName=build-docker-${tag1}-${tag2}


# Build docker image
docker build -t $tag $baseScriptDir

# Run EH build script in docker container
docker run --rm -v $rootDir:$dockerRoot -i $tag bash -<<EOF
set -o errexit
set -o nounset

cd $dockerRoot
mkdir -p $buildDirName
cd $buildDirName
cmake ..
make -j8
EOF

echo
echo "New binaries in ${rootDir}/${buildDirName}/install/bin"
echo

