file(GLOB SOURCES "*.cpp")
add_library(region_analysis ${SOURCES})
target_link_libraries(region_analysis common input reads classification region_spec genotyping alignment filtering ${htslib_static} ${zlib_static})
add_subdirectory(tests)