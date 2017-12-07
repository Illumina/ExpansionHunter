add_executable(graph_test graph_test.cc)
target_link_libraries(graph_test graphs gtest_main)
add_test(NAME graph_test COMMAND graph_test)