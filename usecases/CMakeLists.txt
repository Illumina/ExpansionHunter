file(GLOB SOURCES "*.cc")
add_library(reads_operations ${SOURCES})
target_link_libraries(reads_operations reads)
add_subdirectory(unit_tests)