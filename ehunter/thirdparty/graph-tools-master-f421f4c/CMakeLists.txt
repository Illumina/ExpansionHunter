cmake_minimum_required(VERSION 3.5.0)

project(graphtools CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(BUILD_TESTS OFF CACHE BOOL "Should unit tests be built")

set(USE_ASAN OFF CACHE BOOL "Use clang address sanitizer")
set(USE_MSAN OFF CACHE BOOL "Use clang memory sanitizer")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG")

if (CMAKE_SYSTEM_PROCESSOR MATCHES "^x86_64$")
    if (NOT GRAPHTOOLS_AVX2)
        set(GRAPHTOOLS_VECTORIZATION "-msse4.2")
    else (NOT GRAPHTOOLS_AVX2)
        set(GRAPHTOOLS_VECTORIZATION "-mavx2")
    endif (NOT GRAPHTOOLS_AVX2)
endif (CMAKE_SYSTEM_PROCESSOR MATCHES "^x86_64$")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(GRAPHTOOLS_CXX_OPTIMIZATION_FLAGS "${GRAPHTOOLS_VECTORIZATION} -O2 -ftree-vectorize -finline-functions -fpredictive-commoning -fgcse-after-reload -funswitch-loops -ftree-slp-vectorize -fvect-cost-model -fipa-cp-clone -ftree-phiprop")
endif ()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GRAPHTOOLS_CXX_OPTIMIZATION_FLAGS} -Wall -Werror -pedantic -Wsign-compare -Wno-missing-braces -Wno-deprecated-copy -Wno-unknown-warning-option")

if (USE_ASAN)
    SET(CLANG_ASAN "-O1 -g -fsanitize=address -fno-omit-frame-pointer")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLANG_ASAN}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CLANG_ASAN}")
    SET(BUILD_TESTS ON)
endif ()

if (USE_MSAN)
    SET(CLANG_MSAN "-O1 -g -fsanitize=memory -fno-omit-frame-pointer")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLANG_MSAN}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CLANG_MSAN}")
    SET(BUILD_TESTS ON)
endif ()

#find_package(Threads REQUIRED)
find_package(Boost 1.5 COMPONENTS program_options filesystem system REQUIRED)

file(GLOB SOURCES "src/graphalign/*.cpp" "src/graphalign/*/*.cpp" "src/graphcore/*.cpp" "src/graphutils/*.cpp" "src/graphio/*.cpp")
add_library(graphtools ${SOURCES})
target_include_directories(graphtools PUBLIC "include")
target_include_directories(graphtools PUBLIC "external/include")
target_include_directories(graphtools PUBLIC "external/sparsepp/sparsepp-e40d7a0")
target_link_libraries(graphtools Boost::boost)
target_compile_features(graphtools PRIVATE cxx_range_for)

if (BUILD_TESTS)
    enable_testing()
    # Download and unpack googletest at configure time
    include(GetGoogleTest)
    add_subdirectory(tests)
endif (BUILD_TESTS)
