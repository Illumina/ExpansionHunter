file(GLOB SOURCES "*.cpp")
add_library(stats ${SOURCES})
target_link_libraries(stats common)
add_subdirectory(tests)
