file(GLOB SOURCES "*.cpp")
add_library(region_spec ${SOURCES})
target_link_libraries(region_spec common graphtools)