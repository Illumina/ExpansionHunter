file(GLOB SOURCES "*.cc")
add_library(region_spec ${SOURCES})
target_link_libraries(region_spec common graphtools)
add_subdirectory(unit_tests)