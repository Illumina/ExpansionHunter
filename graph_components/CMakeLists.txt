add_library(graph_components
        ReadClassifier.hh ReadClassifier.cpp GraphVariantAlignmentStats.cpp GraphVariantAlignmentStats.hh)

target_link_libraries(graph_components reads common)

add_subdirectory(tests)