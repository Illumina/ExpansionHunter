file(GLOB sources *.cc)
add_library(common ${sources})
target_link_libraries(common ${Boost_LIBRARIES})