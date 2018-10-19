file(GLOB SOURCES "*.cc")
add_library(reads ${SOURCES})
target_link_libraries(reads graphtools common)
add_subdirectory(unit_tests)