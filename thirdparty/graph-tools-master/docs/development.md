# Building

The following is required to build graph-tools library from source:

- A C++ 11 compliant compiler (tested with GCC 6.3.0),
- CMake 3.5.0 or above.
- BOOST 1.5 or above

To build the library, perform the standard out-of-source CMake build. The unit
tests are not built by default. To build the unit tests, pass `-DBUILD_TESTS=ON`
to CMake.

To also build the included graphIO library, pass `-DBUILD_GRAPHIO=ON` to CMake.
In that case htslib is required. By default htslib is downloaded from github and build during the CMake configure step. Alternatively
set $HTSLIB_INSTALL_PATH to the path (install prefix) of an already installed htslib.
GraphIO also includes a (header-only) version of the nlohmann/json library.

## Incorporating the library into other CMake projects

Copy graph-tools to the source tree of your project. Assuming that the library
was placed in `thirdparty/graph-tools` directory in the project's root, add
`add_subdirectory(thirdparty/graph-tools)` and
`target_link_libraries(<YOUR_PROJECT> graphtools)` to your CMakeLists.txt file.

# Code standards
C++ code should follow the [c++ core guidelines](https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md). 
Code has to pass cppcheck; see [docker-cppcheck.sh](/src/sh/docker-cppcheck.sh) for details.

## Style guide
Naming of identifiers should match the surrounding code. 
Formatting should (exactly) match the output of clang-format 5.0.1; see [format-everything.sh](/src/sh/format-everything.sh) for details.
