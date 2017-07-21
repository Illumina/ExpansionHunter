add_executable(count_test count_test.cc)
target_link_libraries(count_test stats gmock_main gtest)
add_test(NAME count_test COMMAND count_test)