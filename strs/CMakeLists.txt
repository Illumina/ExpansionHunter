add_library(strs
        ReadClassifier.hh ReadClassifier.cpp)

target_link_libraries(strs common reads)

add_subdirectory(tests)