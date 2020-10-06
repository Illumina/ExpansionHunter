file(GLOB SOURCES "*.cpp")
add_library(sample_analysis ${SOURCES})
target_link_libraries(sample_analysis region_analysis common)
add_subdirectory(tests)
