file(GLOB SOURCES "*.cc")
add_library(classification ${SOURCES})
target_link_libraries(classification graphs)
add_subdirectory(unit_tests)