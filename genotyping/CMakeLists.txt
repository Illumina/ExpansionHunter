file(GLOB SOURCES "*.cc")
add_library(genotyping ${SOURCES})
target_link_libraries(genotyping common)
add_subdirectory(unit_tests)
