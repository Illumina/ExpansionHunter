file(GLOB SOURCES "*.cc")
add_library(common ${SOURCES})
target_link_libraries(common ${htslib_static} ${zlib_static} ${Boost_LIBRARIES})
add_subdirectory(unit_tests)
