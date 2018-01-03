file(GLOB SOURCES "*.cc")
add_library(reads ${SOURCES})
target_link_libraries(reads graphs)
add_subdirectory(unit_tests)