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

#include "sample/HtsStreamingReadPairQueue.hh"

namespace ehunter
{

bool HtsStreamingReadPairQueue::insertReadPair(const unsigned locusIndex, ReadPair readPair)
{
    auto& locusAnalyzerQueue(queues_[locusIndex]);
    std::unique_lock<std::mutex> locusAnalyzerQueueLock(locusAnalyzerQueue.mutex);
    const bool wasInActive(not locusAnalyzerQueue.isActive);
    if (wasInActive)
    {
        std::unique_lock<std::mutex> globalLock(mutex_);
        while (activeLocusAnalyzerQueues_ >= maxActiveLocusAnalyzerQueues_)
        {
            cv_.wait(globalLock);
        }
        activeLocusAnalyzerQueues_++;
        globalLock.unlock();

        locusAnalyzerQueue.isActive = true;
    }
    locusAnalyzerQueue.queue.emplace(std::move(readPair));
    return wasInActive;
}

void HtsStreamingReadPairQueue::getNextReadPair(const unsigned locusIndex, boost::optional<ReadPair>& readPair)
{
    auto& locusAnalyzerQueue(queues_[locusIndex]);
    std::unique_lock<std::mutex> locusAnalyzerQueueLock(locusAnalyzerQueue.mutex);

    if (locusAnalyzerQueue.queue.empty())
    {
        std::unique_lock<std::mutex> globalLock(mutex_);
        activeLocusAnalyzerQueues_--;
        globalLock.unlock();

        locusAnalyzerQueue.isActive = false;
        locusAnalyzerQueueLock.unlock();

        cv_.notify_one();
        readPair = boost::none;
    }
    else
    {
        readPair = std::move(locusAnalyzerQueue.queue.front());
        locusAnalyzerQueue.queue.pop();
    }
}

}
