//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Chris Saunders <csaunders@illumina.com>
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

#pragma once

#include <condition_variable>
#include <mutex>
#include <queue>

#include "boost/optional.hpp"

#include "core/Read.hh"
#include "locus/LocusAnalyzer.hh"
#include "sample/AnalyzerFinder.hh"

namespace ehunter
{

/// \brief Custom concurrent queue manager for EH streaming mode.
///
/// THe parallelization strategy used by EH streaming mode has a constraint to have no more than one thread operating on
/// each LocusAnalyzer at a time. This object assists by holding a queue of work items (ReadPairs) for each
/// LocusAnalyzer, managing parallel read/write requests to each queue, and limiting the total number of queues which
/// will be saved.
///
class HtsStreamingReadPairQueue
{
public:
    /// \param[in] maxActiveLocusAnalyzerQueues The max number of non-empty LocusAnalyzer queues to store before
    /// blocking additional input
    ///
    HtsStreamingReadPairQueue(const unsigned maxActiveLocusAnalyzerQueues, const unsigned locusAnalyzerCount)
        : maxActiveLocusAnalyzerQueues_(maxActiveLocusAnalyzerQueues)
        , activeLocusAnalyzerQueues_(0)
        , queues_(locusAnalyzerCount)
    {
    }

    struct ReadPair
    {
        ReadPair(const ReadPair&) = delete;
        ReadPair& operator=(const ReadPair&) = delete;
        ReadPair(ReadPair&& other) = default;
        ReadPair& operator=(ReadPair&& other) = default;

        locus::RegionType regionType;
        AnalyzerInputType inputType;
        Read read;
        Read mate;
    };

    /// \brief Insert a new read pair into the \p locusIndex queue
    ///
    /// If \p locusIndex corresponds to an inactive queue, this will block until the queue can be activated without
    /// exceeding maxActiveLocusAnalyzerQueues.
    ///
    /// \return True if the locusAnalyzer at \p locusIndex was marked as inactive before this method call
    ///
    bool insertReadPair(unsigned locusIndex, ReadPair readPair);

    /// \brief Retrieve a read pair from the \p locusIndex queue
    ///
    /// \param[out] readPair Next read pair enqueued for \p locusIndex, or none if the queue is empty
    ///
    void getNextReadPair(unsigned locusIndex, boost::optional<ReadPair>& readPair);

private:
    struct LocusAnalyzerQueue
    {
        std::queue<ReadPair> queue;

        /// Queue mutex is used to synchronize between the two queue interaction threads (reader and writer)
        std::mutex mutex;

        /// True if a thread is either processing or scheduled to process this queue already
        bool isActive = false;
    };

    const unsigned maxActiveLocusAnalyzerQueues_;
    unsigned activeLocusAnalyzerQueues_;
    std::vector<LocusAnalyzerQueue> queues_;

    /// Global mutex/cv are used to enforce the max active LocusAnalyzerQueue limit
    std::mutex mutex_;
    std::condition_variable cv_;
};

}