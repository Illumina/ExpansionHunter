file(GLOB SOURCES "*.cc")
add_library(stats ${SOURCES})
target_link_libraries(stats common)
add_subdirectory(unit_tests)
