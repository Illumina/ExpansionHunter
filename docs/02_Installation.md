# Installation

Expansion Hunter is designed for Linux and macOS operating systems. A compiled
binary for the latest release can be downloaded from
[here](https://github.com/Illumina/ExpansionHunter/releases). If you wish to
build the program from source follow the instructions below.

## Building from source

Prerequisites:

 - A recent version of [gcc](https://gcc.gnu.org/) or
   [clang](http://clang.llvm.org/) compiler supporting the C++11 standard
   - The minimum gcc version is 5.1
 - [CMake](https://cmake.org/) version 3.13.0 or above
 - Additional development packages, which depend on the operating system:
     - Centos8
       - `bzip2-devel libcurl-devel libstdc++-static openssl-devel xz-devel zlib-devel`
     - Ubuntu 20.04
       - `libbz2-dev libcurl4-openssl-dev liblzma-dev libssl-dev zlib1g-dev `
     - macOS 10.15
       - `xz` (from homebrew)

If the above prerequisites are satisfied, you are ready to
build the program. Note that during the build procedure, cmake will
attempt to download and install `abseil`, `boost`, `googletest`, `htslib`,
and `spdlog` so an active internet connection is required. Assuming
that the source code is contained in a directory `ExpansionHunter/`,
the build procedure can be initiated as follows:

```bash
$ cd ExpansionHunter
$ mkdir build
$ cd build
$ cmake ..
$ make
```

If all the above steps were successful, the ExpansionHunter executable can be found in:

    build/install/bin/ExpansionHunter
