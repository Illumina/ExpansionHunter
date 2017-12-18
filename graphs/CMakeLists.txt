file(GLOB SOURCES "*.cc")
add_library(graphs ${SOURCES})
target_link_libraries(graphs common)
add_subdirectory(unit_tests)