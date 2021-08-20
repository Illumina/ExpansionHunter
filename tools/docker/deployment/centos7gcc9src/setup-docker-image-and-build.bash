#!/usr/bin/env bash

set -o errexit nounset

scriptDir=$(cd $(dirname $0) && pwd -P)
rootDir=$(cd $scriptDir && cd ../../../.. && pwd -P)
dockerRoot=/builder
tag1=ehdeploy
tag2=centos7gcc9src
tag=${tag1}:${tag2}
buildDirName=build-docker-${tag1}-${tag2}


# Build docker image
sudo docker build -t $tag $scriptDir

# Run EH build script in docker container
sudo docker run -v $rootDir:$dockerRoot -i $tag bash -<<EOF
set -o errexit nounset

# build
cd $dockerRoot
mkdir -p $buildDirName
cd $buildDirName
cmake ..
make -j8
EOF

echo
echo "New binaries in ${rootDir}/${buildDirName}/install/bin"
echo

