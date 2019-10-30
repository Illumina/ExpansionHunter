file(GLOB SOURCES "*.cpp")
add_library(alignment ${SOURCES})
target_link_libraries(alignment common graphtools input locus_spec)

add_subdirectory(tests)
