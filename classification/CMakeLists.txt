file(GLOB SOURCES "*.cpp")
add_library(classification ${SOURCES})
target_link_libraries(classification region_spec input graphtools reads)
add_subdirectory(tests)