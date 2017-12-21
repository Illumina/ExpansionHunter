add_executable(read_test read_test.cc)

target_link_libraries(read_test reads gtest_main)

add_test(NAME read_test COMMAND read_test)
