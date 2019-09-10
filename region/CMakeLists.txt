add_library(region STATIC
        Region.hh Region.cpp
        GraphRegion.hh GraphRegion.cpp
        StrFeature.cpp StrFeature.hh
        SmallVariantFeature.cpp SmallVariantFeature.hh
        StrAnalyzer.cpp StrAnalyzer.hh
        VariantAnalyzer.cpp VariantAnalyzer.hh
        VariantFindings.cpp VariantFindings.hh WorkflowBuilder.cpp WorkflowBuilder.hh)

target_link_libraries(region reads alignment filtering region_spec common)