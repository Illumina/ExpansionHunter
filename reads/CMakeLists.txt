file(GLOB SOURCES "*.cpp")
add_library(reads ${SOURCES})
target_link_libraries(reads graphtools common)
add_subdirectory(tests)