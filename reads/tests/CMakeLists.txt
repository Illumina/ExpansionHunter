add_executable(ReadTest ReadTest.cpp)
target_link_libraries(ReadTest reads gtest_main)
add_test(NAME ReadTest COMMAND ReadTest)
