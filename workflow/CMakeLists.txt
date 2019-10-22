add_library(workflow STATIC
        RegionModel.hh RegionModel.cpp
        GraphModel.hh GraphModel.cpp
        GraphStr.cpp GraphStr.hh
        GraphSmallVariant.cpp GraphSmallVariant.hh
        GraphStrAnalyzer.cpp GraphStrAnalyzer.hh
        FeatureAnalyzer.cpp FeatureAnalyzer.hh
        VariantFindings.cpp VariantFindings.hh
        WorkflowBuilder.cpp WorkflowBuilder.hh
        LocusAnalyzer.cpp LocusAnalyzer.hh
        LocusFindings.cpp LocusFindings.hh
        GraphLocusAnalyzer.cpp GraphLocusAnalyzer.hh
        CnvLocusAnalyzer.cpp CnvLocusAnalyzer.hh
        LinearModel.cpp LinearModel.hh
        Feature.cpp Feature.hh
        LinearFeature.cpp LinearFeature.hh
        ReadCountAnalyzer.cpp ReadCountAnalyzer.hh
        GraphVariantAnalyzer.cpp GraphVariantAnalyzer.hh
        IrrPairDetector.cpp IrrPairDetector.hh
        GraphSmallVariantAnalyzer.cpp GraphSmallVariantAnalyzer.hh
        CnvVariantAnalyzer.cpp CnvVariantAnalyzer.hh
        ReadCounter.cpp ReadCounter.hh
        GraphFeature.cpp GraphFeature.hh)

target_link_libraries(workflow strs classification reads alignment filtering region_spec common)