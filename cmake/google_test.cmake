set(GTEST_DIR "../../software/googletest/" CACHE PATH "Google Test path.")

add_subdirectory(${GTEST_DIR} ${CMAKE_BINARY_DIR}/gtest)
include_directories(SYSTEM ${GTEST_DIR}/googlemock/include/ ${GTEST_DIR}/googletest/include/)

function(add_google_test target)
  add_executable(${target} ${ARGN})
  target_link_libraries(${target} gmock_main)
  add_test(${target} ${target})
endfunction()
