file(GLOB SOURCES "*.cc")
add_library(alignment ${SOURCES})
target_link_libraries(alignment common graphtools region_spec)

add_subdirectory(unit_tests)
