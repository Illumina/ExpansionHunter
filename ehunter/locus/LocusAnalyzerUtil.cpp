//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "locus/LocusAnalyzerUtil.hh"

#include <atomic>
#include <thread>

#include "spdlog/spdlog.h"

// Note that ostr.h must be included after spdlog.h
#include "spdlog/fmt/ostr.h"

namespace ehunter
{
namespace locus
{
namespace
{

struct LocusInitThreadSharedData
{
    LocusInitThreadSharedData()
        : isWorkerThreadException(false)
        , locusIndex(0)
    {
    }

    std::atomic<bool> isWorkerThreadException;
    std::atomic<unsigned> locusIndex;
};

/// \brief Data isolated to each locus-initialization thread
///
struct LocusThreadLocalData
{
    std::exception_ptr threadExceptionPtr = nullptr;
};

/// \brief Initialize a series of locus analyzers on one thread
///
void initializeLocusAnalyzerThread(
    const int threadIndex, const RegionCatalog& regionCatalog, const HeuristicParameters& heuristicParams,
    AlignWriterPtr bamletWriter, std::vector<std::unique_ptr<LocusAnalyzer>>& locusAnalyzers,
    LocusInitThreadSharedData& locusInitThreadSharedData,
    std::vector<LocusThreadLocalData>& locusInitThreadLocalDataPool)
{
    LocusThreadLocalData& locusThreadData(locusInitThreadLocalDataPool[threadIndex]);
    std::string locusId = "Unknown";

    try
    {
        const unsigned size(regionCatalog.size());
        while (true)
        {
            if (locusInitThreadSharedData.isWorkerThreadException.load())
            {
                return;
            }
            const auto locusIndex(locusInitThreadSharedData.locusIndex.fetch_add(1));
            if (locusIndex >= size)
            {
                return;
            }

            const auto& locusSpec(regionCatalog[locusIndex]);
            locusId = locusSpec.locusId();

            locusAnalyzers[locusIndex].reset(new LocusAnalyzer(locusSpec, heuristicParams, bamletWriter));
        }
    }
    catch (const std::exception& e)
    {
        locusInitThreadSharedData.isWorkerThreadException.store(true);
        locusThreadData.threadExceptionPtr = std::current_exception();

        spdlog::error(
            "Exception caught in thread {}  while initializing locus: {} : {}", threadIndex, locusId, e.what());
        throw;
    }
    catch (...)
    {
        locusInitThreadSharedData.isWorkerThreadException.store(true);
        locusThreadData.threadExceptionPtr = std::current_exception();

        spdlog::error("Exception caught in thread {}  while initializing locus: {}", threadIndex, locusId);
        throw;
    }
}

}

std::vector<std::unique_ptr<LocusAnalyzer>> initializeLocusAnalyzers(
    const RegionCatalog& regionCatalog, const HeuristicParameters& heuristicParams, AlignWriterPtr bamletWriter,
    const int threadCount)
{
    assert(threadCount >= 1);

    std::vector<std::unique_ptr<LocusAnalyzer>> locusAnalyzers;
    locusAnalyzers.resize(regionCatalog.size());

    LocusInitThreadSharedData locusInitThreadSharedData;
    std::vector<LocusThreadLocalData> locusInitThreadLocalDataPool(threadCount);

    std::vector<std::thread> initThreads;
    for (int threadIndex(0); threadIndex < threadCount; threadIndex++)
    {
        initThreads.emplace_back(
            initializeLocusAnalyzerThread, threadIndex, std::cref(regionCatalog), std::cref(heuristicParams),
            bamletWriter, std::ref(locusAnalyzers), std::ref(locusInitThreadSharedData),
            std::ref(locusInitThreadLocalDataPool));
    }

    // Rethrow exceptions from worker pool in thread order:
    if (locusInitThreadSharedData.isWorkerThreadException.load())
    {
        for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
        {
            const auto& locusThreadData(locusInitThreadLocalDataPool[threadIndex]);
            if (locusThreadData.threadExceptionPtr)
            {
                std::rethrow_exception(locusThreadData.threadExceptionPtr);
            }
        }
    }

    for (auto& initThread : initThreads)
    {
        initThread.join();
    }

    return locusAnalyzers;
}

}
}
