add_library(region STATIC
        RegionModel.hh RegionModel.cpp
        GraphModel.hh GraphModel.cpp
        StrFeature.cpp StrFeature.hh
        SmallVariantFeature.cpp SmallVariantFeature.hh
        StrAnalyzer.cpp StrAnalyzer.hh
        VariantAnalyzer.cpp VariantAnalyzer.hh
        VariantFindings.cpp VariantFindings.hh WorkflowBuilder.cpp WorkflowBuilder.hh LocusAnalyzer.cpp LocusAnalyzer.hh LocusFindings.cpp LocusFindings.hh GraphLocus.cpp GraphLocus.hh CountingModel.cpp CountingModel.hh ModelFeature.cpp ModelFeature.hh GraphFeature.cpp GraphFeature.hh)

target_link_libraries(region reads alignment filtering region_spec common)