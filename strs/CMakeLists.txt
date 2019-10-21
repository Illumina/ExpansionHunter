add_library(strs
        ReadClassifier.hh ReadClassifier.cpp StrAlignmentStats.cpp StrAlignmentStats.hh)

target_link_libraries(strs common reads)

add_subdirectory(tests)