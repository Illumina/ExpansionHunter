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

#include "sample/HtsStreamingSampleAnalysis.hh"

#include <memory>

#include "absl/container/flat_hash_set.h"
#include "spdlog/spdlog.h"
#include <boost/optional.hpp>

#include "core/HtsHelpers.hh"
#include "core/ThreadPool.hh"
#include "locus/LocusAnalyzer.hh"
#include "locus/LocusAnalyzerUtil.hh"
#include "sample/GenomeQueryCollection.hh"
#include "sample/HtsFileStreamer.hh"
#include "sample/HtsStreamingReadPairQueue.hh"

using ehunter::locus::initializeLocusAnalyzers;
using ehunter::locus::LocusAnalyzer;
using graphtools::AlignmentWriter;
using std::string;
using std::vector;

namespace ehunter
{

namespace
{

/// \brief Mutable data shared by all LocusAnalyzer-processing threads
///
class LocusAnalyzerThreadSharedData
{
public:
    LocusAnalyzerThreadSharedData(const unsigned maxActiveLocusAnalyzerQueues, const unsigned locusAnalyzerCount)
        : isWorkerThreadException(false)
        , readPairQueue(maxActiveLocusAnalyzerQueues, locusAnalyzerCount)
    {
    }

    std::atomic<bool> isWorkerThreadException;
    HtsStreamingReadPairQueue readPairQueue;
    vector<std::unique_ptr<LocusAnalyzer>> locusAnalyzers;
};

/// \brief Data isolated to each LocusAnalyzer-processing thread
///
struct LocusAnalyzerThreadLocalData
{
    std::exception_ptr threadExceptionPtr = nullptr;
    std::shared_ptr<graphtools::AlignerSelector> alignerSelectorPtr;
};

/// Process queue for a single LocusAnalyzer on one thread
void processLocusAnalyzerQueue(
    const int threadIndex, LocusAnalyzerThreadSharedData& locusAnalyzerThreadSharedData,
    std::vector<LocusAnalyzerThreadLocalData>& locusAnalyzerThreadLocalDataPool, const unsigned locusIndex)
{
    if (locusAnalyzerThreadSharedData.isWorkerThreadException.load())
    {
        return;
    }

    LocusAnalyzerThreadLocalData& locusAnalyzerThreadData(locusAnalyzerThreadLocalDataPool[threadIndex]);
    auto& locusAnalyzer(*locusAnalyzerThreadSharedData.locusAnalyzers[locusIndex]);

    boost::optional<HtsStreamingReadPairQueue::ReadPair> readPair;

    try
    {
        while (true)
        {
            locusAnalyzerThreadSharedData.readPairQueue.getNextReadPair(locusIndex, readPair);
            if (not readPair)
            {
                break;
            }
            processAnalyzerBundleReadPair(
                locusAnalyzer, readPair->regionType, readPair->inputType, readPair->read, readPair->mate,
                *locusAnalyzerThreadData.alignerSelectorPtr);
        }
    }
    catch (const std::exception& e)
    {
        locusAnalyzerThreadSharedData.isWorkerThreadException.store(true);
        locusAnalyzerThreadData.threadExceptionPtr = std::current_exception();

        std::ostringstream oss;
        oss << "Exception caught in thread " << threadIndex << " while processing read pair queue for locus: `"
            << locusAnalyzer.locusId() << "`";
        if (readPair)
        {
            oss << " current readPair: `" << readPair->read.fragmentId() << "`";
        }
        oss << ": " << e.what();
        spdlog::error(oss.str());
        throw;
    }
    catch (...)
    {
        locusAnalyzerThreadSharedData.isWorkerThreadException.store(true);
        locusAnalyzerThreadData.threadExceptionPtr = std::current_exception();

        std::ostringstream oss;
        oss << "Exception caught in thread " << threadIndex << " while processing read pair queue for locus: `"
            << locusAnalyzer.locusId() << "`";
        if (readPair)
        {
            oss << " current readPair: `" << readPair->read.fragmentId() << "`";
        }
        spdlog::error(oss.str());
        throw;
    }
}

/// \brief Mutable data shared by all SampleFindings-processing threads
///
class SampleFindingsThreadSharedData
{
public:
    SampleFindingsThreadSharedData()
        : isWorkerThreadException(false)
        , locusIndex(0)
    {
    }

    std::atomic<bool> isWorkerThreadException;
    std::atomic<unsigned> locusIndex;
};

/// \brief Data isolated to each SampleFindings-processing thread
///
struct SampleFindingsThreadLocalData
{
    std::exception_ptr threadExceptionPtr = nullptr;
};

/// \brief Analyze a series of loci on one thread
///
void analyzeLocus(
    const int threadIndex, const Sex sampleSex, vector<std::unique_ptr<LocusAnalyzer>>& locusAnalyzers,
    SampleFindings& sampleFindings, SampleFindingsThreadSharedData& sampleFindingsThreadSharedData,
    std::vector<SampleFindingsThreadLocalData>& sampleFindingsThreadLocalData)
{
    SampleFindingsThreadLocalData& sampleFindingsThreadData(sampleFindingsThreadLocalData[threadIndex]);
    std::string locusId = "Unknown";

    try
    {

        const unsigned size(locusAnalyzers.size());
        while (true)
        {
            if (sampleFindingsThreadSharedData.isWorkerThreadException.load())
            {
                return;
            }
            const auto locusIndex(sampleFindingsThreadSharedData.locusIndex.fetch_add(1));
            if (locusIndex >= size)
            {
                return;
            }

            auto& locusAnalyzer(*locusAnalyzers[locusIndex]);
            locusId = locusAnalyzer.locusId();
            sampleFindings[locusIndex] = locusAnalyzer.analyze(sampleSex, boost::none);
        }
    }
    catch (const std::exception& e)
    {
        sampleFindingsThreadSharedData.isWorkerThreadException = true;
        sampleFindingsThreadData.threadExceptionPtr = std::current_exception();

        spdlog::error("Exception caught in thread {} while analyzing locus: {} : {}", threadIndex, locusId, e.what());
        throw;
    }
    catch (...)
    {
        sampleFindingsThreadSharedData.isWorkerThreadException = true;
        sampleFindingsThreadData.threadExceptionPtr = std::current_exception();

        spdlog::error("Unknown exception caught in thread {} while analyzing locus: {}", threadIndex, locusId);
        throw;
    }
}

}

SampleFindings htsStreamingSampleAnalysis(
    const InputPaths& inputPaths, Sex sampleSex, const HeuristicParameters& heuristicParams, const int threadCount,
    const RegionCatalog& regionCatalog, locus::AlignWriterPtr bamletWriter)
{
    // Setup thread-specific data structures and thread pool
    const unsigned maxActiveLocusAnalyzerQueues(threadCount + 5);
    const unsigned locusAnalyzerCount(regionCatalog.size());
    LocusAnalyzerThreadSharedData locusAnalyzerThreadSharedData(maxActiveLocusAnalyzerQueues, locusAnalyzerCount);
    std::vector<LocusAnalyzerThreadLocalData> locusAnalyzerThreadLocalDataPool(threadCount);
    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        auto& locusAnalyzerThreadData(locusAnalyzerThreadLocalDataPool[threadIndex]);
        locusAnalyzerThreadData.alignerSelectorPtr.reset(
            new graphtools::AlignerSelector(heuristicParams.alignerType()));
    }
    ctpl::thread_pool pool(threadCount);

    spdlog::info("Initializing all loci");
    graphtools::AlignerSelector alignerSelector(heuristicParams.alignerType());
    locusAnalyzerThreadSharedData.locusAnalyzers
        = initializeLocusAnalyzers(regionCatalog, heuristicParams, bamletWriter, threadCount);
    GenomeQueryCollection genomeQuery(locusAnalyzerThreadSharedData.locusAnalyzers);

    spdlog::info("Streaming reads");

    auto ReadHash = [](const Read& read) { return std::hash<std::string>()(read.fragmentId()); };
    auto ReadEq = [](const Read& read1, const Read& read2) { return (read1.fragmentId() == read2.fragmentId()); };
    using ReadCatalog = absl::flat_hash_set<Read, decltype(ReadHash), decltype(ReadEq)>;
    ReadCatalog unpairedReads(1000, ReadHash, ReadEq);

    const unsigned htsDecompressionThreads(std::min(threadCount, 12));
    htshelpers::HtsFileStreamer readStreamer(inputPaths.htsFile(), inputPaths.reference(), htsDecompressionThreads);
    while (readStreamer.trySeekingToNextPrimaryAlignment() && readStreamer.isStreamingAlignedReads())
    {
        // Stop processing reads if an exception is thrown in the worker pool:
        if (locusAnalyzerThreadSharedData.isWorkerThreadException.load())
        {
            break;
        }

        const bool isReadNearTargetRegion = genomeQuery.targetRegionMask.query(
            readStreamer.currentReadContigId(), readStreamer.currentReadPosition());
        const bool isMateNearTargetRegion = genomeQuery.targetRegionMask.query(
            readStreamer.currentMateContigId(), readStreamer.currentMatePosition());
        if (!isReadNearTargetRegion && !isMateNearTargetRegion)
        {
            continue;
        }

        if (not readStreamer.currentIsPaired())
        {
            continue;
        }

        Read read = readStreamer.decodeRead();
        const auto mateIterator = unpairedReads.find(read);
        if (mateIterator == unpairedReads.end())
        {
            unpairedReads.emplace(std::move(read));
            continue;
        }
        Read mate = std::move(*mateIterator);
        unpairedReads.erase(mateIterator);

        const int64_t readEnd = readStreamer.currentReadPosition() + read.sequence().length();
        const int64_t mateEnd = readStreamer.currentMatePosition() + mate.sequence().length();

        vector<AnalyzerBundle> analyzerBundles = genomeQuery.analyzerFinder.query(
            readStreamer.currentReadContigId(), readStreamer.currentReadPosition(), readEnd,
            readStreamer.currentMateContigId(), readStreamer.currentMatePosition(), mateEnd);

        const unsigned bundleCount(analyzerBundles.size());
        for (unsigned bundleIndex(0); bundleIndex < bundleCount; ++bundleIndex)
        {
            auto& bundle(analyzerBundles[bundleIndex]);
            auto sendReadPair = [&](HtsStreamingReadPairQueue::ReadPair readPair)
            {
                if (locusAnalyzerThreadSharedData.readPairQueue.insertReadPair(bundle.locusIndex, std::move(readPair)))
                {
                    pool.push(
                        processLocusAnalyzerQueue, std::ref(locusAnalyzerThreadSharedData),
                        std::ref(locusAnalyzerThreadLocalDataPool), bundle.locusIndex);
                }
            };

            if ((bundleIndex + 1) < bundleCount)
            {
                sendReadPair({ bundle.regionType, bundle.inputType, read, mate });
            }
            else
            {
                sendReadPair({ bundle.regionType, bundle.inputType, std::move(read), std::move(mate) });
            }
        }
    }

    pool.stop(true);

    // Rethrow exceptions from the pool in thread order:
    if (locusAnalyzerThreadSharedData.isWorkerThreadException.load())
    {
        for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
        {
            const LocusAnalyzerThreadLocalData& locusAnalyzerThreadData(locusAnalyzerThreadLocalDataPool[threadIndex]);
            if (locusAnalyzerThreadData.threadExceptionPtr)
            {
                std::rethrow_exception(locusAnalyzerThreadData.threadExceptionPtr);
            }
        }
    }

    spdlog::info("Analyzing read evidence");

    SampleFindingsThreadSharedData sampleFindingsThreadSharedData;
    std::vector<SampleFindingsThreadLocalData> sampleFindingsThreadLocalDataPool(threadCount);

    const unsigned locusCount(locusAnalyzerThreadSharedData.locusAnalyzers.size());
    SampleFindings sampleFindings(locusCount);

    // Start all sampleFindings worker threads
    std::vector<std::thread> sampleFindingsThreads;
    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        sampleFindingsThreads.emplace_back(
            analyzeLocus, threadIndex, sampleSex, std::ref(locusAnalyzerThreadSharedData.locusAnalyzers),
            std::ref(sampleFindings), std::ref(sampleFindingsThreadSharedData),
            std::ref(sampleFindingsThreadLocalDataPool));
    }

    // Rethrow exceptions from worker pool in thread order:
    if (sampleFindingsThreadSharedData.isWorkerThreadException.load())
    {
        for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
        {
            const auto& sampleFindingsThreadData(sampleFindingsThreadLocalDataPool[threadIndex]);
            if (sampleFindingsThreadData.threadExceptionPtr)
            {
                std::rethrow_exception(sampleFindingsThreadData.threadExceptionPtr);
            }
        }
    }

    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        sampleFindingsThreads[threadIndex].join();
    }

    return sampleFindings;
}

}
