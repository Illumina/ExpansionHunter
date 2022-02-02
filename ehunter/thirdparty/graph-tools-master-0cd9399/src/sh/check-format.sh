#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SOURCE_DIR="${DIR}/../../src"
export CLANG_FORMAT=clang-format-5.0

set +e
DIFFS=0
for x in $(find ${SOURCE_DIR} -iname *.h -o -iname *.hh -o -iname *.cpp); do
    ${CLANG_FORMAT} -style=file ${x} | diff --ignore-blank-lines --strip-trailing-cr ${x} -
    if [[ $? != 0 ]]; then
        echo "Differences found in $x"
        DIFFS=1
    fi
done

if [[ ${DIFFS} == 0 ]]; then
    echo "No formatting issues found."
    exit 0
else
    echo "Some formatting issues were found."
    exit 1
fi
