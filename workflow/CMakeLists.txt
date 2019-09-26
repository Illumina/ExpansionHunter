add_library(workflow STATIC
        RegionModel.hh RegionModel.cpp
        GraphModel.hh GraphModel.cpp
        GraphStr.cpp GraphStr.hh
        GraphSmallVariant.cpp GraphSmallVariant.hh
        GraphStrAnalyzer.cpp GraphStrAnalyzer.hh
        FeatureAnalyzer.cpp FeatureAnalyzer.hh
        VariantFindings.cpp VariantFindings.hh WorkflowBuilder.cpp WorkflowBuilder.hh LocusAnalyzer.cpp LocusAnalyzer.hh LocusFindings.cpp LocusFindings.hh GraphLocusAnalyzer.cpp GraphLocusAnalyzer.hh LinearModel.cpp LinearModel.hh RegionModelFeature.cpp RegionModelFeature.hh GraphVariant.cpp GraphVariant.hh LinearModelFeature.cpp LinearModelFeature.hh StatsAnalyzer.cpp StatsAnalyzer.hh GraphVariantAnalyzer.cpp GraphVariantAnalyzer.hh OfftargetFeature.cpp OfftargetFeature.hh GraphSmallVariantAnalyzer.cpp GraphSmallVariantAnalyzer.hh)

target_link_libraries(workflow strs reads alignment filtering region_spec common)