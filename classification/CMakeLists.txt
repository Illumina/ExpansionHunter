file(GLOB SOURCES "*.cc")
add_library(classification ${SOURCES})
target_link_libraries(classification region_spec graphtools reads)
add_subdirectory(unit_tests)