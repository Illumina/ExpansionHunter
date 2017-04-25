file(GLOB SOURCES "*.cc")
add_library(genotyping ${SOURCES})
add_subdirectory(unit_tests)
