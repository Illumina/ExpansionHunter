file(GLOB SOURCES "*.cc")
add_library(sample_analysis ${SOURCES})
target_link_libraries(sample_analysis region_analysis common)
