#!/usr/bin/env bash

set -o errexit nounset

dockerRoot=/builder
scriptDir=$(cd $(dirname $0) && pwd -P)
rootDir=$(cd $scriptDir && cd ../../../.. && pwd -P)
buildDirName=build-docker
tag=ehdebug:mallinfo2


# Build docker image
sudo docker build -t $tag .

# Run EH build script in docker container
sudo docker run -v $rootDir:$dockerRoot -i $tag bash -<<EOF 
set -o errexit nounset

cd $dockerRoot 
mkdir -p $buildDirName 
cd $buildDirName 
cmake .. 
make -j8
EOF

echo
echo "New binaries in ${rootDir}/${buildDirName}/install/bin"
echo

