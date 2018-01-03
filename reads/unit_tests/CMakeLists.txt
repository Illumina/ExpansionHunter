add_executable(read_test read_test.cc)
add_executable(read_operations_test read_operations_test.cc)

target_link_libraries(read_test reads gtest_main)
target_link_libraries(read_operations_test reads gtest_main)

add_test(NAME read_test COMMAND read_test)
add_test(NAME read_operations_test COMMAND read_operations_test)
