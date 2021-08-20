######################### Google Test ############################
# Download and unpack googletest at configure time

FILE(WRITE "${CMAKE_BINARY_DIR}/external/googletest-build/CMakeLists.txt" "\
cmake_minimum_required(VERSION 2.8.12)
project(googletest-build NONE)
include(ExternalProject)
ExternalProject_Add(googletest
        # GIT_REPOSITORY    https://github.com/google/googletest.git
        # GIT_TAG           release-1.8.0
        URL               \"${CMAKE_SOURCE_DIR}/external/googletest-release-1.10.0.tar.gz\"
        URL_HASH          MD5=ecd1fa65e7de707cd5c00bdac56022cd
        SOURCE_DIR        \"${CMAKE_BINARY_DIR}/external/googletest-src\"
        BINARY_DIR        \"${CMAKE_BINARY_DIR}/external/googletest-build\"
        CONFIGURE_COMMAND \"\"
        BUILD_COMMAND     \"\"
        INSTALL_COMMAND   \"\"
        TEST_COMMAND      \"\"
        )"
)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/googletest-build )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/googletest-build )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/external/googletest-src
        ${CMAKE_BINARY_DIR}/external/googletest-build)
#add_subdirectory(${CMAKE_BINARY_DIR}/external/googletest-src/googletest
#        ${CMAKE_BINARY_DIR}/external/googletest-build-test)
#add_subdirectory(${CMAKE_BINARY_DIR}/external/googletest-src/googlemock
#        ${CMAKE_BINARY_DIR}/external/googletest-build-mock)

##################################################################