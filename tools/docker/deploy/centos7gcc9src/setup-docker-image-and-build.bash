#!/usr/bin/env bash

scriptDir=$(cd "$(dirname "$0")" && pwd -P)

${scriptDir}/../../shared/setup-docker-image-and-build.bash "5.8" "$0"

