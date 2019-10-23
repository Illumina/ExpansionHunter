
add_library(classification
        AlignmentSummary.hh AlignmentSummary.cpp
        SmallVariantAlignmentClassifier.hh SmallVariantAlignmentClassifier.cpp
        StrAlignmentClassifier.hh StrAlignmentClassifier.cpp
        )
target_link_libraries(classification reads locus_spec input alignment graphtools)
add_subdirectory(tests)