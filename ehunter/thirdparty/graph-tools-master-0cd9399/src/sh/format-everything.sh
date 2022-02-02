#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CLANG_FORMAT=clang-format

if [ ! -x "$(command -v ${CLANG_FORMAT})" ]; then
    CLANG_FORMAT=clang-format-5.0
fi

if [ ! -x "$(command -v ${CLANG_FORMAT})" ]; then
    CLANG_FORMAT=clang-format-mp-5.0
fi

if [ ! -x "$(command -v ${CLANG_FORMAT})" ]; then
    echo "Clang-format not found, please install as clang-format or clang-format-mp-5.0"
    exit 1
fi

if [[ -z "$(${CLANG_FORMAT} --version | grep 'version 5')" ]]; then
    echo "Clang-format has the wrong version, please install as clang-format version 5.0"
    exit 1
fi

find ${DIR}/../../src -iname *.h -o -iname *.hh -o -iname *.cpp | xargs ${CLANG_FORMAT} -i -style=file
find ${DIR}/../../include -iname *.h -o -iname *.hh -o -iname *.cpp | xargs ${CLANG_FORMAT} -i -style=file
find ${DIR}/../../tests -iname *.h -o -iname *.hh -o -iname *.cpp | xargs ${CLANG_FORMAT} -i -style=file
