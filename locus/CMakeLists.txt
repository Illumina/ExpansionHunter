add_library(locus STATIC
        Locus.hh Locus.cpp
        GraphLocus.hh GraphLocus.cpp StrFeature.cpp StrFeature.hh SmallVariantFeature.cpp SmallVariantFeature.hh StrAnalyzer.cpp StrAnalyzer.hh VariantAnalyzer.cpp VariantAnalyzer.hh VariantFindings.cpp VariantFindings.hh)

target_link_libraries(locus reads alignment filtering common)