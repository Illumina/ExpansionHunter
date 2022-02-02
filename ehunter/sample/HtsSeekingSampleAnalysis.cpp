//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Chris Saunders <csaunders@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//

#include "sample/HtsSeekingSampleAnalysis.hh"

#include <atomic>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>
#include <boost/smart_ptr/make_unique.hpp>

// clang-format off
// Note that spdlog.h must be included before ostr.h
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
// clang-format on

#include "core/ReadPairs.hh"
#include "locus/LocusAnalyzer.hh"
#include "sample/AnalyzerFinder.hh"
#include "sample/HtsFileSeeker.hh"
#include "sample/IndexBasedDepthEstimate.hh"
#include "sample/MateExtractor.hh"

using boost::make_unique;
using boost::optional;
using ehunter::htshelpers::HtsFileSeeker;
using ehunter::locus::LocusAnalyzer;
using graphtools::AlignmentWriter;
using std::ostream;
using std::string;
using std::unique_ptr;
using std::unordered_map;
using std::vector;

namespace ehunter
{

namespace
{
using AlignmentStatsCatalog = unordered_map<ReadId, LinearAlignmentStats, boost::hash<ReadId>>;

vector<GenomicRegion>
combineRegions(const vector<GenomicRegion>& targetRegions, const vector<GenomicRegion>& offtargetRegions)
{
    vector<GenomicRegion> combinedRegions(targetRegions);
    combinedRegions.insert(combinedRegions.end(), offtargetRegions.begin(), offtargetRegions.end());
    return combinedRegions;
}

bool checkIfMatesWereMappedNearby(const LinearAlignmentStats& alignmentStats)
{
    const int kMaxMateDistance = 1000;
    if ((alignmentStats.chromId == alignmentStats.mateChromId)
        && (std::abs(alignmentStats.pos - alignmentStats.matePos) < kMaxMateDistance))
    {
        return true;
    }
    return false;
}

void recoverMates(
    htshelpers::MateExtractor& mateExtractor, AlignmentStatsCatalog& alignmentStatsCatalog, ReadPairs& readPairs)
{
    for (auto& fragmentIdAndReadPair : readPairs)
    {
        ReadPair& readPair = fragmentIdAndReadPair.second;

        if (readPair.numMatesSet() == 2)
        {
            continue;
        }

        const Read& read = readPair.firstMate ? *readPair.firstMate : *readPair.secondMate;

        const auto alignmentStatsIterator = alignmentStatsCatalog.find(read.readId());
        if (alignmentStatsIterator == alignmentStatsCatalog.end())
        {
            throw std::logic_error("Cannot recover mate of uncatalogued read");
        }
        const LinearAlignmentStats& alignmentStats = alignmentStatsIterator->second;

        if (!checkIfMatesWereMappedNearby(alignmentStats))
        {
            LinearAlignmentStats mateStats;
            optional<Read> optionalMate = mateExtractor.extractMate(read, alignmentStats, mateStats);
            if (optionalMate)
            {
                const Read& mate = *optionalMate;
                alignmentStatsCatalog.emplace(std::make_pair(mate.readId(), alignmentStats));
                readPairs.AddMateToExistingRead(mate);
            }
            else
            {
                spdlog::warn("Could not recover the mate of {}", read.readId());
            }
        }
    }
}

ReadPairs collectCandidateReads(
    const vector<GenomicRegion>& targetRegions, const vector<GenomicRegion>& offtargetRegions,
    AlignmentStatsCatalog& alignmentStatsCatalog, HtsFileSeeker& htsFileSeeker,
    htshelpers::MateExtractor& mateExtractor)
{
    vector<GenomicRegion> regionsWithReads = combineRegions(targetRegions, offtargetRegions);
    ReadPairs readPairs;

    for (const auto& regionWithReads : regionsWithReads)
    {
        const int numReadsBeforeCollection = readPairs.NumReads();
        htsFileSeeker.setRegion(regionWithReads);
        while (htsFileSeeker.trySeekingToNextPrimaryAlignment())
        {
            LinearAlignmentStats alignmentStats;
            Read read = htsFileSeeker.decodeRead(alignmentStats);
            if (alignmentStats.isPaired)
            {
                alignmentStatsCatalog.emplace(std::make_pair(read.readId(), alignmentStats));
                readPairs.Add(std::move(read));
            }
            else
            {
                spdlog::warn("Skipping {} because it is unpaired", read.readId());
            }
        }
        const int numReadsCollected = readPairs.NumReads() - numReadsBeforeCollection;
        spdlog::debug("Collected {} reads from {}", numReadsCollected, regionWithReads);
    }

    const int numReadsBeforeRecovery = readPairs.NumReads();
    recoverMates(mateExtractor, alignmentStatsCatalog, readPairs);
    const int numReadsAfterRecovery = readPairs.NumReads() - numReadsBeforeRecovery;
    spdlog::debug("Recovered {} reads", numReadsAfterRecovery);

    return readPairs;
}

void analyzeReadPair(
    vector<unique_ptr<LocusAnalyzer>>& locusAnalyzers, AnalyzerFinder& analyzerFinder, Read& read, Read& mate,
    const AlignmentStatsCatalog& alignmentStats, graphtools::AlignerSelector& alignerSelector)
{
    const auto readStatsIter = alignmentStats.find(read.readId());
    const auto mateStatsIter = alignmentStats.find(mate.readId());

    if (readStatsIter == alignmentStats.end() || mateStatsIter == alignmentStats.end())
    {
        throw std::logic_error("Could not to find alignment stats for " + read.fragmentId());
    }

    const LinearAlignmentStats& readStats = readStatsIter->second;
    const LinearAlignmentStats& mateStats = mateStatsIter->second;

    const int64_t readEnd = readStats.pos + read.sequence().length();
    const int64_t mateEnd = mateStats.pos + mate.sequence().length();
    vector<AnalyzerBundle> analyzers
        = analyzerFinder.query(readStats.chromId, readStats.pos, readEnd, mateStats.chromId, mateStats.pos, mateEnd);

    if (analyzers.empty())
    {
        return;
    }

    assert(analyzers.size() == 1);
    auto& analyzer(analyzers.front());
    processAnalyzerBundleReadPair(
        *locusAnalyzers[analyzer.locusIndex], analyzer.regionType, analyzer.inputType, read, mate, alignerSelector);
}

void analyzeRead(
    vector<unique_ptr<LocusAnalyzer>>& locusAnalyzers, AnalyzerFinder& analyzerFinder, Read& read,
    const AlignmentStatsCatalog& alignmentStats, graphtools::AlignerSelector& alignerSelector)
{
    const auto readStatsIter = alignmentStats.find(read.readId());

    if (readStatsIter == alignmentStats.end())
    {
        throw std::logic_error("Could not to find alignment stats for " + read.fragmentId());
    }

    const LinearAlignmentStats& readStats = readStatsIter->second;
    const int64_t readEnd = readStats.pos + read.sequence().length();

    vector<AnalyzerBundle> analyzers = analyzerFinder.query(readStats.chromId, readStats.pos, readEnd);

    if (analyzers.empty())
    {
        return;
    }

    assert(analyzers.size() == 1);
    auto& analyzer(analyzers.front());
    locusAnalyzers[analyzer.locusIndex]->processMates(read, nullptr, analyzer.regionType, alignerSelector);
}

void processReads(
    vector<unique_ptr<LocusAnalyzer>>& locusAnalyzers, ReadPairs& candidateReadPairs,
    const AlignmentStatsCatalog& alignmentStats, AnalyzerFinder& analyzerFinder,
    graphtools::AlignerSelector& alignerSelector)
{
    for (auto& fragmentIdAndReads : candidateReadPairs)
    {
        auto& readPair = fragmentIdAndReads.second;
        if (readPair.numMatesSet() == 2)
        {
            analyzeReadPair(
                locusAnalyzers, analyzerFinder, *readPair.firstMate, *readPair.secondMate, alignmentStats,
                alignerSelector);
        }
        else
        {
            Read& read = readPair.firstMate ? *readPair.firstMate : *readPair.secondMate;
            analyzeRead(locusAnalyzers, analyzerFinder, read, alignmentStats, alignerSelector);
        }
    }
}

/// \brief Mutable data shared by all worker threads
///
class LocusThreadSharedData
{
public:
    LocusThreadSharedData()
        : isWorkerThreadException(false)
        , locusIndex(0)
    {
    }

    std::atomic<bool> isWorkerThreadException;
    std::atomic<unsigned> locusIndex;
};

/// \brief Data isolated to each locus-processing thread
///
struct LocusThreadLocalData
{
    std::exception_ptr threadExceptionPtr = nullptr;
};

/// \brief Process a series of loci on one thread
///
void processLocus(
    const int threadIndex, const InputPaths& inputPaths, const Sex sampleSex,
    const HeuristicParameters& heuristicParams, const RegionCatalog& regionCatalog,
    locus::AlignWriterPtr alignmentWriter, SampleFindings& sampleFindings, LocusThreadSharedData& locusThreadSharedData,
    std::vector<LocusThreadLocalData>& locusThreadLocalDataPool)
{
    LocusThreadLocalData& locusThreadData(locusThreadLocalDataPool[threadIndex]);
    std::string locusId = "Unknown";

    try
    {
        HtsFileSeeker htsFileSeeker(inputPaths.htsFile(), inputPaths.reference());
        htshelpers::MateExtractor mateExtractor(inputPaths.htsFile(), inputPaths.reference());
        graphtools::AlignerSelector alignerSelector(heuristicParams.alignerType());

        const unsigned size(regionCatalog.size());
        while (true)
        {
            if (locusThreadSharedData.isWorkerThreadException.load())
            {
                return;
            }
            const auto locusIndex(locusThreadSharedData.locusIndex.fetch_add(1));
            if (locusIndex >= size)
            {
                return;
            }

            const auto& locusSpec(regionCatalog[locusIndex]);
            locusId = locusSpec.locusId();

            spdlog::info("Analyzing {}", locusId);
            vector<unique_ptr<LocusAnalyzer>> locusAnalyzers;
            auto analyzer(make_unique<LocusAnalyzer>(locusSpec, heuristicParams, alignmentWriter));
            locusAnalyzers.emplace_back(std::move(analyzer));
            AnalyzerFinder analyzerFinder(locusAnalyzers);

            AlignmentStatsCatalog alignmentStats;
            ReadPairs readPairs = collectCandidateReads(
                locusSpec.targetReadExtractionRegions(), locusSpec.offtargetReadExtractionRegions(), alignmentStats,
                htsFileSeeker, mateExtractor);

            processReads(locusAnalyzers, readPairs, alignmentStats, analyzerFinder, alignerSelector);

            sampleFindings[locusIndex] = locusAnalyzers.front()->analyze(sampleSex, boost::none);
        }
    }
    catch (const std::exception& e)
    {
        locusThreadSharedData.isWorkerThreadException = true;
        locusThreadData.threadExceptionPtr = std::current_exception();

        spdlog::error("Exception caught in thread {} while processing locus: {} : {}", threadIndex, locusId, e.what());
        throw;
    }
    catch (...)
    {
        locusThreadSharedData.isWorkerThreadException = true;
        locusThreadData.threadExceptionPtr = std::current_exception();

        spdlog::error("Unknown exception caught in thread {} while processing locus: {}", threadIndex, locusId);
        throw;
    }
}
}

SampleFindings htsSeekingSampleAnalysis(
    const InputPaths& inputPaths, Sex sampleSex, const HeuristicParameters& heuristicParams, const int threadCount,
    const RegionCatalog& regionCatalog, locus::AlignWriterPtr alignmentWriter)
{
    if (ehunter::isURL(inputPaths.htsFile()))
    {
        // For URL input paths, the index needs to be downloaded in advance if seeking mode is using multiple threads.
        // This is needed because htslib has no protection against the race condition created by multiple threads
        // independently downloading this index to the same file path.
        //
        if (threadCount > 1)
        {
            (void)HtsFileSeeker(inputPaths.htsFile(), inputPaths.reference());
        }
    }

    LocusThreadSharedData locusThreadSharedData;
    std::vector<LocusThreadLocalData> locusThreadLocalDataPool(threadCount);

    const unsigned locusCount(regionCatalog.size());
    SampleFindings sampleFindings(locusCount);

    // Start all locus worker threads
    std::vector<std::thread> locusThreads;
    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        locusThreads.emplace_back(
            processLocus, threadIndex, std::cref(inputPaths), sampleSex, std::cref(heuristicParams),
            std::cref(regionCatalog), alignmentWriter, std::ref(sampleFindings), std::ref(locusThreadSharedData),
            std::ref(locusThreadLocalDataPool));
    }

    // Rethrow exceptions from worker pool in thread order:
    if (locusThreadSharedData.isWorkerThreadException.load())
    {
        for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
        {
            const auto& locusThreadData(locusThreadLocalDataPool[threadIndex]);
            if (locusThreadData.threadExceptionPtr)
            {
                std::rethrow_exception(locusThreadData.threadExceptionPtr);
            }
        }
    }

    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        locusThreads[threadIndex].join();
    }

    return sampleFindings;
}

}
