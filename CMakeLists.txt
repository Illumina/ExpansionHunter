cmake_minimum_required(VERSION 3.1.0)
project(ExpansionHunter CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

enable_testing()

include(ExternalProject)

option(BUILD_TESTS "Controls if unit tests are build" OFF)
if(BUILD_TESTS)
  include(google_test)
endif()

ExternalProject_Add(zlib
	PREFIX ${CMAKE_BINARY_DIR}/thirdparty/zlib
	GIT_REPOSITORY "https://github.com/madler/zlib.git"
	GIT_TAG "v1.2.8"
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
	CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/thirdparty/zlib/src/zlib/configure --prefix=${CMAKE_BINARY_DIR}/thirdparty/zlib --static
	INSTALL_DIR ${CMAKE_BINARY_DIR}/thirdparty/zlib
	LOG_DOWNLOAD 1
	LOG_INSTALL 1
)

ExternalProject_Add(htslib
	PREFIX ${CMAKE_BINARY_DIR}/thirdparty/htslib
        GIT_REPOSITORY "https://github.com/samtools/htslib.git"
	GIT_TAG "1.3.1"
        UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
	BUILD_COMMAND make
	INSTALL_COMMAND make install prefix=${CMAKE_BINARY_DIR}/thirdparty/htslib
	LOG_DOWNLOAD 1
)

include_directories(${CMAKE_BINARY_DIR}/thirdparty/zlib/include)
set(zlib_static ${CMAKE_BINARY_DIR}/thirdparty/zlib/lib/libz.a)
set(htslib_static ${CMAKE_BINARY_DIR}/thirdparty/htslib/lib/libhts.a)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost 1.4 COMPONENTS program_options filesystem regex date_time system REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${CMAKE_BINARY_DIR}/thirdparty/htslib/include)

add_subdirectory(purity)
add_subdirectory(rep_align)

file(GLOB SOURCES "src/*.cc")

add_executable(ExpansionHunter ${SOURCES})
target_compile_features(ExpansionHunter PRIVATE cxx_range_for)

add_dependencies(htslib zlib)
add_dependencies(ExpansionHunter htslib)

target_link_libraries(ExpansionHunter purity rep_align pthread ${htslib_static} ${zlib_static} ${Boost_LIBRARIES})
