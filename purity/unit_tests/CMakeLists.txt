add_executable(purity_test purity_test.cc)
target_link_libraries(purity_test purity gtest_main)
add_test(NAME purity_test COMMAND purity_test)