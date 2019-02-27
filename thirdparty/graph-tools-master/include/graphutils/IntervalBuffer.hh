//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Peter Krusche <pkrusche@illumina.com>
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

#pragma once

#include <cstddef>
#include <cstdint>
#include <list>
#include <memory>
#include <utility>

namespace intervals
{

class IntervalBuffer
{
public:
    /** tracks intervals over a number of lanes */
    IntervalBuffer();
    IntervalBuffer(IntervalBuffer const& rhs);
    virtual ~IntervalBuffer();
    IntervalBuffer& operator=(IntervalBuffer const& rhs);

    /**
     * @brief Add an interval to a lane
     *
     * @param start interval coordinates
     * @param end interval coordinates
     * @param lane lane to add to
     */
    void addInterval(int64_t start, int64_t end, size_t lane);

    /**
     * @brief Advance buffer, discarding all intervals with end < to
     *
     * @param to interval minimum end position; pass -1 to clear buffer
     */
    void advance(int64_t to);

    /**
     * @brief Check if interval is fully covered in a given lane
     */
    bool isCovered(int64_t start, int64_t end, size_t lane) const;

    /**
     * @brief Check if interval is partially covered in a given lane
     */
    bool hasOverlap(int64_t start, int64_t end, size_t lane) const;

    /**
     * Get the intervals for a particular lane
     * @return intervals for a particular lane
     */
    std::list<std::pair<int, int>> getIntervals(size_t lane) const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};
}
