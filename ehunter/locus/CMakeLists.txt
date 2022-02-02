target_sources(ExpansionHunterLib # Requires CMake 3.13 or later
        PRIVATE
        AlignmentBuffer.hh AlignmentBuffer.cpp
        IrrPairFinder.hh IrrPairFinder.cpp
        LocusAligner.hh LocusAligner.cpp
        LocusAnalyzer.hh LocusAnalyzer.cpp
        LocusAnalyzerUtil.hh LocusAnalyzerUtil.cpp
        LocusFindings.hh LocusFindings.cpp
        LocusSpecification.hh LocusSpecification.cpp
        RepeatAnalyzer.hh RepeatAnalyzer.cpp
        RFC1MotifAnalysis.hh RFC1MotifAnalysis.cpp
        RFC1MotifAnalysisUtil.hh RFC1MotifAnalysisUtil.cpp
        RFC1Status.hh
        SmallVariantAnalyzer.hh SmallVariantAnalyzer.cpp
        VariantAnalyzer.hh VariantAnalyzer.cpp
        VariantFindings.hh VariantFindings.cpp
        VariantSpecification.hh VariantSpecification.cpp
        )

target_sources(UnitTests # Requires CMake 3.13 or later
        PRIVATE
        IrrPairFinderTest.cpp
        LocusAlignerTest.cpp
        LocusAnalyzerTest.cpp
        )
