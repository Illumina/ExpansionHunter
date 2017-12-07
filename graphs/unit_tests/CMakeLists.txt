add_executable(graph_test graph_test.cc)
add_executable(graph_builders_test graph_builders_test.cc)

target_link_libraries(graph_test graphs gtest_main)
target_link_libraries(graph_builders_test graphs gtest_main)

add_test(NAME graph_test COMMAND graph_test)
add_test(NAME graph_builders_test COMMAND graph_builders_test)
