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

/**
 * \brief Store a list of non-intersecting intervals
 *
 *
 * \file IntervalList.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <list>
#include <map>
#include <vector>

namespace intervals
{

/** interval interface */
struct interval
{
    explicit interval(int64_t _start = -1, int64_t _end = -1)
        : start(_start)
        , end(_end)
    {
    }
    virtual ~interval() = default;

    // merge two intervals
    virtual void merge(interval const& rhs)
    {
        if (start < 0)
        {
            start = rhs.start;
        }
        else
        {
            start = std::min(rhs.start, start);
        }
        if (end < 0)
        {
            end = rhs.end;
        }
        else
        {
            end = std::max(rhs.end, end);
        }
    }

    virtual void resize(int64_t _start, int64_t _end)
    {
        if (_start >= 0)
        {
            start = _start;
        }
        if (_end >= 0)
        {
            end = _end;
        }
    }

    int64_t start, end;
};

template <class interval_t = interval> class IntervalList
{
public:
    virtual ~IntervalList() = default;

    void add(interval_t const& iv)
    {
        if (iv.start > iv.end)
        {
            return;
        }

        // find first interval that ends after start - 1 (so we join adjacent intervals)
        auto ivr = intervals.lower_bound(iv.start - 1);
        if (ivr == intervals.end())
        {
            // no interval ends after start => prepend
            intervals[iv.end] = iv;
        }
        else
        {
            // first interval that ends after start
            // check overlap

            // do they overlap
            // we know ivr->first >= start
            if (ivr->second.start <= iv.end)
            {
                // overlap => merge
                interval_t tmp = ivr->second;
                tmp.merge(iv);
                intervals.erase(ivr);
                add(tmp);
            }
            else
            {
                // no overlap:
                // x   y        s       f
                //-[---]--------[-------]-------
                //       |----|
                //    start  end
                // (y must be < start because otherwise we'd have found it with lower_bound)
                // => insert interval
                intervals[iv.end] = iv;
            }
        }
    }

    /** get all intervals that overlap a range */
    interval_t query(int64_t start, int64_t end) const
    {
        interval_t ivs;

        if (end < start)
        {
            return ivs;
        }

        // find first interval that ends after start
        auto it = intervals.cbegin();

        int x = 0;
        while (it != intervals.cend() && it->second.end < start && x < 3)
        {
            ++it;
            ++x;
        }

        // if using advance, we don't need to search -- the first interval will
        // already be the one we're looking for
        if (it != intervals.cend() && it->second.end < start)
        {
            it = intervals.lower_bound(start);
        }

        // overlap if it->second.start <= end
        while (it != intervals.cend() &&
               // check if overlap
               it->second.start <= end)
        {
            ivs.merge(it->second);
            ++it;
        }
        return ivs;
    }

    /** get all intervals that overlap a range */
    template <class container_t = std::list<interval_t>> void get(int64_t start, int64_t end, container_t& output)
    {
        if (end < start)
        {
            return;
        }

        // find first interval that ends after start
        auto it = intervals.lower_bound(start);

        // overlap if it->second.start <= end
        while (it != intervals.end() &&
               // check if overlap
               it->second.start <= end)
        {
            output.push_back(it->second);
            ++it;
        }
    }

    /** reset / remove interval range */
    void remove_from(int64_t start)
    {
        // find first interval that ends after start
        auto it = intervals.lower_bound(start);

        // keep all intervals that end before start
        // know:
        // it->second.start <= it->second.end &&
        //            start <= it->second.end &&
        //            start <= end
        // region to delete starts inside interval?
        // (if interval starts exactly at start, we can remove it)
        if (it != intervals.end() && it->second.start < start)
        {
            interval_t tmp = it->second;
            intervals.erase(it, intervals.end());
            // start is > end if end == -1
            tmp.resize(tmp.start, start - 1);
            intervals[tmp.end] = tmp;
        }
        else
        {
            // remove stuff that ends afterwards, if any
            intervals.erase(it, intervals.end());
        }
    }

    /** reset / remove interval range */
    void remove_to(int64_t end)
    {
        if (end < 0)
        {
            intervals.clear();
            return;
        }

        // find first interval that ends after end
        auto it = intervals.lower_bound(end);
        if (it != intervals.end() && it->second.start <= end)
        {
            if (it->first > end)
            {
                it->second.resize(end + 1, it->first);
            }
            else // fully contained => also remove
            {
                ++it;
            }
        }
        intervals.erase(intervals.begin(), it);
    }

    /** remove values outside region */
    void keep_only(int64_t _start, int64_t _end)
    {
        remove_to(_start - 1);
        remove_from(_end + 1);
    }

    /** return intervals in arbitrary container */
    template <class container_t = std::list<interval_t>> container_t getIntervals() const
    {
        container_t t;
        for (auto const& iv : intervals)
        {
            t.push_back(iv.second);
        }
        return t;
    }

protected:
    std::map<int64_t, interval_t> intervals;

    void ensure_is_interval()
    {
        // make sure interval_t is derived from interval
        interval* iv = new interval_t();
        delete iv;
    }
};
}
