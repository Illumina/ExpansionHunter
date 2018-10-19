file(GLOB SOURCES "*.cc")
add_library(region_analysis ${SOURCES})
target_link_libraries(region_analysis common reads classification region_spec genotyping alignment filtering ${htslib_static} ${zlib_static})
add_subdirectory(unit_tests)