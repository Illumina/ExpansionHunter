# Installation

Expansion Hunter is designed for Linux and macOS operating systems.
A compiled binary for the latest release can be downloaded from 
[here](https://github.com/Illumina/ExpansionHunter/releases). If 
you wish to build the program from source follow the instructions 
below.

Building from source
--------------------

Prerequisites:

 - A recent version of [GCC](https://gcc.gnu.org/) or 
   [clang](http://clang.llvm.org/) compiler supporting C++11.
 - [CMake](https://cmake.org/) version 3.2.0 or above.
 - [Boost C++ Libraries](http://www.boost.org/) version 1.57.0 or
   above.

If you the above prerequisites are satisfied, you are ready to
build the program. Note that during the build procedure, cmake will 
attempt to download and install [HTSlib](http://www.htslib.org) and 
[zlib](https://github.com/madler/zlib) so an active internet 
connection is required. Assuming that the source code is contained 
in a directory `ExpansionHunter/`, the build procedure can be 
initiated as follows:

```bash
$ cd ExpansionHunter
$ mkdir build
$ cd build
$ cmake ..
$ make
```

Note that if Boost is installed in a non-default location then its path should be specified with `BOOST_ROOT` in the cmake command above:

```bash
$ cmake -DBOOST_ROOT=/path/to/boost/ ..
```

If all of the above steps were successful, the `build` directory now contains ExpansionHunter executable.