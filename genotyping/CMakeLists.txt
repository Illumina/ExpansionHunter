file(GLOB SOURCES "*.cpp")
add_library(genotyping ${SOURCES})
target_link_libraries(genotyping common)
add_subdirectory(tests)
