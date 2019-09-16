cmake_minimum_required(VERSION 3.1.0)
project(ExpansionHunter CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

enable_testing()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(ExternalProject)

######################### Google Test ############################
# Download and unpack googletest at configure time
configure_file(cmake/google_test.cmake googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
if (result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif ()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
if (result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif ()

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
        ${CMAKE_BINARY_DIR}/googletest-build)
##################################################################


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

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.4 REQUIRED COMPONENTS program_options filesystem regex date_time system)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(SYSTEM ${Boost_INCLUDE_DIR})
include_directories(${CMAKE_BINARY_DIR}/thirdparty/htslib/include)

add_subdirectory(thirdparty/graph-tools-master)

add_compile_options(-Werror -pedantic -Wall -Wextra)

add_subdirectory(common)
add_subdirectory(genotyping)
add_subdirectory(reads)
add_subdirectory(classification)
add_subdirectory(region_spec)
add_subdirectory(sample_analysis)
add_subdirectory(input)
add_subdirectory(output)
add_subdirectory(alignment)
add_subdirectory(stats)
add_subdirectory(filtering)
add_subdirectory(region)
add_subdirectory(strs)

file(GLOB SOURCES "src/*.cpp")
add_executable(ExpansionHunter ${SOURCES})
target_compile_features(ExpansionHunter PRIVATE cxx_range_for)
target_link_libraries(ExpansionHunter graphtools common region genotyping region_spec sample_analysis input output
        alignment filtering stats classification strs pthread ${Boost_LIBRARIES})
install(TARGETS ExpansionHunter DESTINATION bin)

add_dependencies(htslib zlib)
add_dependencies(common htslib)
