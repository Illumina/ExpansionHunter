file(GLOB SOURCES "*.cc")
add_library(graphs ${SOURCES})
add_subdirectory(unit_tests)