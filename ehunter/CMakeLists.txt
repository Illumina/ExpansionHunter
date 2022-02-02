cmake_minimum_required(VERSION 3.13)
project(ExpansionHunter LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE)
    set(DEFAULT_CMAKE_BUILD_TYPE Release)
    set(CMAKE_BUILD_TYPE ${DEFAULT_CMAKE_BUILD_TYPE} CACHE STRING
            "Choose the type of build (default: ${DEFAULT_CMAKE_BUILD_TYPE})" FORCE)
endif ()

message(STATUS "BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.73 REQUIRED COMPONENTS program_options filesystem system)


add_compile_options(-Wall -Werror -Wextra)

enable_testing()

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL REQUIRED)
find_package(GTest REQUIRED)
find_package(Threads REQUIRED)
find_package(spdlog REQUIRED)
find_package(absl REQUIRED)

if (NOT APPLE)
    find_package(OpenSSL REQUIRED)
endif ()

find_library(htslib libhts.a)
find_library(htslib hts)

add_subdirectory(thirdparty/graph-tools-master-0cd9399)

add_library(ExpansionHunterLib
        alignment/AlignmentClassifier.hh alignment/AlignmentClassifier.cpp
        alignment/AlignmentFilters.hh alignment/AlignmentFilters.cpp
        alignment/ClassifierOfAlignmentsToVariant.hh alignment/ClassifierOfAlignmentsToVariant.cpp
        alignment/GraphVariantAlignmentStats.hh alignment/GraphVariantAlignmentStats.cpp
        alignment/GreedyAlignmentIntersector.hh alignment/GreedyAlignmentIntersector.cpp
        alignment/HighQualityBaseRunFinder.hh alignment/HighQualityBaseRunFinder.cpp
        alignment/OperationsOnAlignments.hh alignment/OperationsOnAlignments.cpp
        alignment/OrientationPredictor.hh alignment/OrientationPredictor.cpp
        alignment/SoftclippingAligner.hh alignment/SoftclippingAligner.cpp
        core/Common.hh core/Common.cpp
        core/ConcurrentQueue.hh
        core/CountTable.hh core/CountTable.cpp
        core/GenomicRegion.hh core/GenomicRegion.cpp
        core/HtsHelpers.hh core/HtsHelpers.cpp
        core/Parameters.hh core/Parameters.cpp
        core/Reference.hh core/Reference.cpp
        core/ReferenceContigInfo.hh core/ReferenceContigInfo.cpp
        core/LocusStats.hh core/LocusStats.cpp
        core/LogSum.hh
        core/Read.hh core/Read.cpp
        core/ReadPairs.hh core/ReadPairs.cpp
        core/ReadSupportCalculator.hh core/ReadSupportCalculator.cpp
        core/ThreadPool.hh
        core/WeightedPurityCalculator.hh core/WeightedPurityCalculator.cpp
        genotyping/AlignMatrix.hh genotyping/AlignMatrix.cpp
        genotyping/AlignMatrixFiltering.hh genotyping/AlignMatrixFiltering.cpp
        genotyping/AlleleChecker.hh genotyping/AlleleChecker.cpp
        genotyping/FragLogliks.hh genotyping/FragLogliks.cpp
        genotyping/OneAlleleStrGenotyper.hh genotyping/OneAlleleStrGenotyper.cpp
        genotyping/RepeatGenotype.hh genotyping/RepeatGenotype.cpp
        genotyping/SmallVariantGenotype.hh genotyping/SmallVariantGenotype.cpp
        genotyping/SmallVariantGenotyper.hh genotyping/SmallVariantGenotyper.cpp
        genotyping/StrAlign.hh genotyping/StrAlign.cpp
        genotyping/StrGenotyper.hh genotyping/StrGenotyper.cpp
        genotyping/TwoAlleleStrGenotyper.hh genotyping/TwoAlleleStrGenotyper.cpp
        io/BamletWriter.hh io/BamletWriter.cpp
        io/CatalogLoading.hh io/CatalogLoading.cpp
        io/GraphBlueprint.hh io/GraphBlueprint.cpp
        io/JsonWriter.hh io/JsonWriter.cpp
        io/LocusSpecDecoding.hh io/LocusSpecDecoding.cpp
        io/ParameterLoading.hh io/ParameterLoading.cpp
        io/RegionGraph.hh io/RegionGraph.cpp
        io/SampleStats.hh io/SampleStats.cpp
        io/VcfHeader.hh io/VcfHeader.cpp
        io/VcfWriter.hh io/VcfWriter.cpp
        io/VcfWriterHelpers.hh io/VcfWriterHelpers.cpp
        sample/AnalyzerFinder.hh sample/AnalyzerFinder.cpp
        sample/GenomeMask.hh sample/GenomeMask.cpp
        sample/GenomeQueryCollection.hh sample/GenomeQueryCollection.cpp
        sample/HtsFileSeeker.hh sample/HtsFileSeeker.cpp
        sample/HtsFileStreamer.hh sample/HtsFileStreamer.cpp
        sample/HtsSeekingSampleAnalysis.hh sample/HtsSeekingSampleAnalysis.cpp
        sample/HtsStreamingReadPairQueue.hh sample/HtsStreamingReadPairQueue.cpp
        sample/HtsStreamingSampleAnalysis.hh sample/HtsStreamingSampleAnalysis.cpp
        sample/IndexBasedDepthEstimate.hh sample/IndexBasedDepthEstimate.cpp
        sample/MateExtractor.hh sample/MateExtractor.cpp
        )


set(CTPL_INCLUDE_DIRS thirdparty/ctpl/ctpl-0.0.2)
target_include_directories(ExpansionHunterLib PUBLIC ${CMAKE_SOURCE_DIR})
target_include_directories(ExpansionHunterLib SYSTEM PUBLIC
        ${Boost_INCLUDE_DIRS}
        ${LIBLZMA_INCLUDE_DIRS}
        ${CURL_INCLUDE_DIRS}
        ${CTPL_INCLUDE_DIRS}
        )

# Set static linking of gcc standard libraries to simplify binary distribution
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(STATIC_FLAGS -static-libgcc -static-libstdc++)
endif ()

target_link_libraries(ExpansionHunterLib PUBLIC
        ${STATIC_FLAGS}
        graphtools
        ${htslib}
        ${Boost_LIBRARIES}
        ${LIBLZMA_LIBRARIES}
        ${CURL_LIBRARIES}
        ZLIB::ZLIB
        BZip2::BZip2
        Threads::Threads
        spdlog::spdlog
        absl::flat_hash_set)

if (NOT APPLE)
    target_link_libraries(ExpansionHunterLib PUBLIC
            OpenSSL::Crypto)
endif ()

add_executable(ExpansionHunter
        app/ExpansionHunter.cpp
        )
target_link_libraries(ExpansionHunter ExpansionHunterLib)

add_executable(UnitTests
        tests/AlignMatrixTest.cpp
        tests/AlignmentClassifierTest.cpp
        tests/AlignmentSummaryTest.cpp
        tests/AlleleCheckerTest.cpp
        tests/ClassifierOfAlignmentsToVariantTest.cpp
        tests/ConcurrentQueueTest.cpp
        tests/CountTableTest.cpp
        tests/FragLogliksTest.cpp
        tests/GenomeMaskTest.cpp
        tests/GenomicRegionTest.cpp
        tests/GraphAlignmentOperationsTest.cpp
        tests/GraphBlueprintTest.cpp
        tests/GreedyAlignmentIntersectorTest.cpp
        tests/HighQualityBaseRunFinderTest.cpp
        tests/LocusStatsTest.cpp
        tests/ReadSupportCalculatorTest.cpp
        tests/ReadTest.cpp
        tests/RegionGraphTest.cpp
        tests/RepeatAnalyzerTest.cpp
        tests/RepeatGenotypeTest.cpp
        tests/RFC1MotifAnalysisUtilTest.cpp
        tests/SmallVariantGenotyperTest.cpp
        tests/SoftclippingAlignerTest.cpp
        tests/StrAlignTest.cpp
        tests/StrGenotyperTest.cpp
        tests/UnitTests.cpp
        tests/WeightedPurityCalculatorTest.cpp
        )
add_subdirectory(locus)

target_link_libraries(UnitTests ExpansionHunterLib GTest::GTest GTest::Main)

add_test(NAME UnitTests COMMAND UnitTests)

install(TARGETS ExpansionHunter UnitTests RUNTIME DESTINATION bin)
