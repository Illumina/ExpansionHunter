file(GLOB SOURCES "*.cpp")
add_library(input ${SOURCES})
target_link_libraries(input common stats region_spec)
add_subdirectory(tests)