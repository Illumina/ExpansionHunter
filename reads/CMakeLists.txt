add_library(reads Read.hh Read.cpp ReadPairs.hh ReadPairs.cpp)
target_link_libraries(reads graphtools common)
add_subdirectory(tests)