cmake_minimum_required(VERSION 3.13)
project(ExpansionHunter)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check for required packages:
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL REQUIRED)

if(NOT APPLE)
  find_package(OpenSSL REQUIRED)
endif()

include(ExternalProject)
set(installDir ${CMAKE_CURRENT_BINARY_DIR}/install)


ExternalProject_Add(htslib
	BUILD_IN_SOURCE YES
	URL https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
	CONFIGURE_COMMAND ./configure --prefix=${installDir} --enable-libcurl --enable-s3
	BUILD_COMMAND $(MAKE) lib-static
)


# Setup boost user config so that it uses the same c++ compiler as other components:
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set (BOOST_COMPILER_TAG "gcc")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	set (BOOST_COMPILER_TAG "clang")
	set(B2_OPTIONS ${B2_OPTIONS} "toolset=clang")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
	set (BOOST_COMPILER_TAG "darwin")
	set(B2_OPTIONS ${B2_OPTIONS} "toolset=clang")
endif()

if (BOOST_COMPILER_TAG)
	set(BOOST_UCONFIG "${CMAKE_BINARY_DIR}/user-config.jam")
	file(WRITE "${BOOST_UCONFIG}" "using ${BOOST_COMPILER_TAG} : : \"${CMAKE_CXX_COMPILER}\" ;\n")
	set(BOOST_PATCH_COMMAND ${CMAKE_COMMAND} -E copy ${BOOST_UCONFIG} tools/build/src/user-config.jam)
endif()

ExternalProject_Add(Boost
	BUILD_IN_SOURCE YES
	URL https://archives.boost.io/release/1.73.0/source/boost_1_73_0.tar.bz2
	PATCH_COMMAND ${BOOST_PATCH_COMMAND}
	CONFIGURE_COMMAND ./bootstrap.sh --prefix=${installDir}
	BUILD_COMMAND ./b2 install -j8  --with-filesystem --with-system --with-program_options link=static ${B2_OPTIONS}
	INSTALL_COMMAND ""
)


ExternalProject_Add(spdlog
	URL https://github.com/gabime/spdlog/archive/refs/tags/v1.6.1.tar.gz
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
)


ExternalProject_Add(googletest
	URL https://github.com/google/googletest/archive/refs/tags/release-1.10.0.tar.gz
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
)


# CMAKE_CXX_STANDARD is required for macOS 10.15
ExternalProject_Add(abseil
	URL https://github.com/abseil/abseil-cpp/archive/refs/tags/20210324.2.tar.gz
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DCMAKE_CXX_STANDARD=11
		-DBUILD_TESTING=OFF
)


ExternalProject_Add(ehunter
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/ehunter
	BUILD_ALWAYS YES
	TEST_COMMAND "ctest"
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_PREFIX_PATH:PATH=${installDir}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${CMAKE_SOURCE_DIR}/ehunter
)


ExternalProject_Add_StepDependencies(ehunter configure
	Boost spdlog htslib googletest abseil)

