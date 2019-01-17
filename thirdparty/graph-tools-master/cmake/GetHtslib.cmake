# Download and build htslib at configure time

set(HTSLIB_INSTALL_PATH "$ENV{HTSLIB_INSTALL_PATH}" CACHE STRING "Specify a pre-built version of htslib (install prefix).")

set(HTSLIB_INCLUDED_VERSION "1.3.1")

if (IS_DIRECTORY ${HTSLIB_INSTALL_PATH})
    message( "Using pre-built htslib from ${HTSLIB_INSTALL_PATH}")
    message( WARNING "htslib <= 1.7 is known to leak memory with cram files")
    set(HTSLIB_PATH ${HTSLIB_INSTALL_PATH})
else()
    message( "Using included htslib" )
    set(HTSLIB_PATH ${CMAKE_BINARY_DIR}/external/htslib-install)

    FILE(WRITE "${CMAKE_BINARY_DIR}/external/htslib-build/CMakeLists.txt" "
        cmake_minimum_required(VERSION 2.8.5)\n
        project(htslib-build NONE)
        include(ExternalProject)
        ExternalProject_Add(htslib
        GIT_REPOSITORY \"https://github.com/samtools/htslib.git\"
        GIT_CONFIG \"http.sslVerify=false\"
        GIT_TAG \"${HTSLIB_INCLUDED_VERSION}\"
        SOURCE_DIR        \"${CMAKE_BINARY_DIR}/external/htslib-src\"
        INSTALL_DIR       \"${HTSLIB_PATH}\"
        CONFIGURE_COMMAND \"\"
        PATCH_COMMAND \"\"
        BUILD_COMMAND make -C <SOURCE_DIR> prefix=<INSTALL_DIR>
        INSTALL_COMMAND make -C <SOURCE_DIR> install prefix=<INSTALL_DIR> )")

    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/htslib-build )
    if(result)
        message(FATAL_ERROR "CMake step for htslib failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/htslib-build )
    if(result)
        message(FATAL_ERROR "Build step for htslib failed: ${result}")
    endif()
endif()

## Define HTSlib library for downstream use

include(FindZLIB)
include(FindBZip2)
include(FindLibLZMA)
add_library(htslib STATIC IMPORTED GLOBAL)
set_property(TARGET htslib PROPERTY INTERFACE_INCLUDE_DIRECTORIES
        ${ZLIB_INCLUDE_DIR} ${BZIP2_INCLUDE_DIR} ${LIBLZMA_INCLUDE_DIRS})
set_property(TARGET htslib APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
        ${HTSLIB_PATH}/include)
set_property(TARGET htslib PROPERTY IMPORTED_LOCATION
        "${HTSLIB_PATH}/lib/libhts.a")
set_property(TARGET htslib PROPERTY INTERFACE_LINK_LIBRARIES
        ${BZIP2_LIBRARIES} ${ZLIB_LIBRARIES})

