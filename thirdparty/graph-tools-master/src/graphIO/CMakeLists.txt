cmake_minimum_required(VERSION 3.5.0)

include(GetHtslib)

set(SOURCES "GraphJson.cpp" "ReferenceGenome.cpp" BamWriter.cpp)
add_library(graphIO ${SOURCES})
target_include_directories(graphIO PUBLIC "../../include")
target_include_directories(graphIO PUBLIC "../../external/include")

target_link_libraries(graphIO Boost::filesystem htslib)