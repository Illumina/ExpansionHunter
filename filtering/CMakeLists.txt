file(GLOB SOURCES "*.cc")
add_library(filtering ${SOURCES})
target_link_libraries(filtering common graphtools)

add_subdirectory(unit_tests)
