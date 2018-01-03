file(GLOB SOURCES "*.cc")
add_library(region_analysis ${SOURCES})
target_link_libraries(region_analysis reads)
add_subdirectory(unit_tests)