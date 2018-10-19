// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2017 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/**
 * \brief Interval buffer implementation
 *
 * \file IntervalBuffer.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "graphutils/IntervalBuffer.hh"
#include "graphutils/IntervalList.hh"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <vector>

namespace intervals
{

struct IntervalBuffer::Impl
{
    Impl() = default;
    Impl(Impl const& rhs) = default;
    Impl& operator=(Impl const& rhs) = default;
    Impl(Impl&& rhs) = delete;
    Impl& operator=(Impl&& rhs) = delete;
    ~Impl() = default;

    typedef IntervalList<interval> ivmap_t;

    std::vector<ivmap_t> lanes;
};

/** tracks intervals over a number of lanes */
IntervalBuffer::IntervalBuffer()
    : pimpl_(new Impl())
{
}

IntervalBuffer::~IntervalBuffer() = default;

IntervalBuffer::IntervalBuffer(IntervalBuffer const& rhs)
    : pimpl_(new Impl(*rhs.pimpl_))
{
}

IntervalBuffer& IntervalBuffer::operator=(IntervalBuffer const& rhs)
{
    if (&rhs == this)
    {
        return *this;
    }
    pimpl_.reset(new Impl(*rhs.pimpl_));
    return *this;
}

/**
 * @brief Add an interval to a lane
 *
 * @param start interval coordinates
 * @param end interval coordinates
 * @param lane lane to add to
 * @return Interval identifier
 */
void IntervalBuffer::addInterval(int64_t start, int64_t end, size_t lane)
{
    if (start > end)
    {
        return;
    }
    if (pimpl_->lanes.size() <= lane)
    {
        pimpl_->lanes.resize(lane + 1);
    }
    pimpl_->lanes[lane].add(interval(start, end));
}

/**
 * @brief Advance buffer, discarding all intervals with start < to
 *
 * @param to interval minimum end position
 */
void IntervalBuffer::advance(int64_t to)
{
    if (to < 0)
    {
        pimpl_->lanes.clear();
        return;
    }

    for (auto& lane : pimpl_->lanes)
    {
        lane.remove_to(to - 1);
    }
}

/**
 * @brief Check if interval is fully covered in a given lane
 */
bool IntervalBuffer::isCovered(int64_t start, int64_t end, size_t lane) const
{
    if (lane >= pimpl_->lanes.size())
    {
        return false;
    }

    // intervals of zero length count as covered
    if (end < start)
    {
        return true;
    }

    std::list<interval> is;
    pimpl_->lanes[lane].get(start, end, is);
    if (is.size() != 1)
    {
        // if we overlap with more than one interval, then there must be a gap
        return false;
    }

    interval& it = is.front();
    return it.start <= start && it.end >= end;
}

/**
 * @brief Check if interval is partially covered in a given lane
 */
bool IntervalBuffer::hasOverlap(int64_t start, int64_t end, size_t lane) const
{
    if (lane >= pimpl_->lanes.size())
    {
        return false;
    }

    // intervals of zero length count as covered
    if (end < start)
    {
        return true;
    }

    interval it = pimpl_->lanes[lane].query(start, end);
    return it.start >= 0 && it.end >= 0 && it.end - it.start + 1 > 0;
}

/**
 * Get the intervals for a particular lane
 * @return intervals for a particular lane
 */
std::list<std::pair<int, int>> IntervalBuffer::getIntervals(size_t lane) const
{
    if (lane >= pimpl_->lanes.size())
    {
        throw std::logic_error(std::string("Unknown lane: ") + std::to_string(lane));
    }
    std::list<std::pair<int, int>> result;
    for (auto const& iv : pimpl_->lanes[lane].getIntervals())
    {
        result.emplace_back(iv.start, iv.end);
    }
    return result;
}
}
