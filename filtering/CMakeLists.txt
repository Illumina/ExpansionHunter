file(GLOB SOURCES "*.cpp")
add_library(filtering ${SOURCES})
target_link_libraries(filtering common graphtools)

add_subdirectory(tests)
