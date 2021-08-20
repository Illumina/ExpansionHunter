cmake_minimum_required(VERSION 3.13)
project(ExpansionHunter)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check for required packages:
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL REQUIRED)
find_package(Git REQUIRED)

include(ExternalProject)
set(installDir ${CMAKE_CURRENT_BINARY_DIR}/install)


ExternalProject_Add(htslib
	BUILD_IN_SOURCE YES
	GIT_REPOSITORY "https://github.com/samtools/htslib.git"
	GIT_TAG "1.10.2"
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ""
	BUILD_COMMAND $(MAKE)
	INSTALL_COMMAND $(MAKE) install prefix=${installDir}
	LOG_DOWNLOAD YES
)


# Setup user config to force same c++ compiler in boost
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
    set(BOOST_UCONFIG "${CMAKE_BINARY_DIR}/user-config.jam")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        file(WRITE "${BOOST_UCONFIG}" "using gcc : : \"${CMAKE_CXX_COMPILER}\" ;\n")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        file(WRITE "${BOOST_UCONFIG}" "using clang : : \"${CMAKE_CXX_COMPILER}\" ;\n")
        set(B2_OPTIONS ${B2_OPTIONS} "toolset=clang")
    endif()
    set(BOOST_PATCH_COMMAND ${CMAKE_COMMAND} -E copy ${BOOST_UCONFIG} tools/build/src/user-config.jam)
endif()

ExternalProject_Add(Boost
	BUILD_IN_SOURCE YES
	URL https://boostorg.jfrog.io/artifactory/main/release/1.73.0/source/boost_1_73_0.tar.bz2
	UPDATE_COMMAND ""
	PATCH_COMMAND ${BOOST_PATCH_COMMAND}
	CONFIGURE_COMMAND ./bootstrap.sh --prefix=${installDir}/lib
	BUILD_COMMAND ./b2 install -j8   --prefix=${installDir} --with-filesystem --with-system --with-program_options link=static ${B2_OPTIONS}
	INSTALL_COMMAND ""
)


ExternalProject_Add(spdlog
	GIT_REPOSITORY "https://github.com/gabime/spdlog.git"
	GIT_TAG "v1.6.1"
	UPDATE_COMMAND ""
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
)


ExternalProject_Add(googletest
	GIT_REPOSITORY "https://github.com/google/googletest.git"
	GIT_TAG "release-1.10.0"
	UPDATE_COMMAND ""
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
)


# CMAKE_CXX_STANDARD is required for OS X 10.15
ExternalProject_Add(abseil
	GIT_REPOSITORY "https://github.com/abseil/abseil-cpp"
	GIT_TAG "20210324.2"
	UPDATE_COMMAND ""
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
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${CMAKE_SOURCE_DIR}/ehunter
)


ExternalProject_Add_StepDependencies(ehunter configure
	Boost spdlog htslib googletest abseil)

