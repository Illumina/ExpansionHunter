file(GLOB SOURCES "*.cpp")
add_library(common ${SOURCES})
target_link_libraries(common region_spec ${htslib_static} ${zlib_static} ${Boost_LIBRARIES})
add_subdirectory(tests)
