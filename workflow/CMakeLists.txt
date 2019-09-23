add_library(workflow STATIC
        RegionModel.hh RegionModel.cpp
        GraphModel.hh GraphModel.cpp
        StrFeature.cpp StrFeature.hh
        SmallVariantFeature.cpp SmallVariantFeature.hh
        StrAnalyzer.cpp StrAnalyzer.hh
        FeatureAnalyzer.cpp FeatureAnalyzer.hh
        VariantFindings.cpp VariantFindings.hh WorkflowBuilder.cpp WorkflowBuilder.hh LocusAnalyzer.cpp LocusAnalyzer.hh LocusFindings.cpp LocusFindings.hh GraphLocusAnalyzer.cpp GraphLocusAnalyzer.hh CountingModel.cpp CountingModel.hh ModelFeature.cpp ModelFeature.hh GraphFeature.cpp GraphFeature.hh CountingFeature.cpp CountingFeature.hh StatsAnalyzer.cpp StatsAnalyzer.hh GraphVariantAnalyzer.cpp GraphVariantAnalyzer.hh PairedIrrFeature.cpp PairedIrrFeature.hh)

target_link_libraries(workflow strs reads alignment filtering region_spec common)