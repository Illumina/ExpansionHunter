file(GLOB SOURCES "*.cc")
add_library(output ${SOURCES})
target_link_libraries(output common stats region_spec)
