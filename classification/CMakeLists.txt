file(GLOB SOURCES "*.cpp")
add_library(classification ${SOURCES})
target_link_libraries(classification reads region_spec input alignment graphtools)
add_subdirectory(tests)