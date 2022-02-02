#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

set -e

PULL=1
while test $# -gt 0
do
    case "$1" in
        --no-pull) echo "Not pulling image."; PULL=0
            ;;
        *) echo "Unknown argument $1"; exit 1
            ;;
    esac
    shift
done

if [[ ${PULL} -ne 0 ]]; then
    docker pull ilmncgrpmi/cpp-dev:latest
fi

LOGFILE=$(mktemp)
docker run -v ${DIR}/../..:/graphtools-source --rm ilmncgrpmi/cpp-dev:latest /bin/bash -c "cd /graphtools-source && \
   cppcheck --enable=all --suppress=missingIncludeSystem --suppress=unusedFunction --inline-suppr -I/graphtools-source/include /graphtools-source/src > /dev/null" &> ${LOGFILE}
   
cat ${LOGFILE}

MESSAGES=$(cat ${LOGFILE} | wc -l)
if [[ ${MESSAGES} -ne 0 ]]; then
    echo "$MESSAGES cppcheck errors were found!"
    exit 1
fi

exit 0
