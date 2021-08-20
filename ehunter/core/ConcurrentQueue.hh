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

namespace ehunter
{

template <typename T> class ConcurrentQueue
{
public:
    void push(T const& data)
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            queue_.push(data);
        }
        cv_.notify_one();
    }

    bool empty() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.empty();
    }

    void pop(T& value)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        while (queue_.empty())
        {
            cv_.wait(lock);
        }

        value = queue_.front();
        queue_.pop();
    }

private:
    std::queue<T> queue_;
    mutable std::mutex mutex_;
    std::condition_variable cv_;
};

}
