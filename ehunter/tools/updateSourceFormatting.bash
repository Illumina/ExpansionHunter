#!/usr/bin/env bash
#
# Run clang-format on all project cxx source using clang-format v7.0.1
#

set -o nounset

this_dir=$(dirname $0)

cformat=$this_dir/clang-format

cxx_base_dir=$(cd $this_dir/..; pwd -P)

get_cpp_files() {
    for sdir in alignment app core genotyping io locus sample tests; do
        find $cxx_base_dir/$sdir -type f \( -name *.cpp -or -name *.hh \)
    done
}

# remove windows line endings from source:
get_cpp_files | xargs -P8 -n1 sed $'s/\r$//' -i

# general c++ source reformatting:
get_cpp_files | xargs -P8 -n1 $cformat -style=file -i
