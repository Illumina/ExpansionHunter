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

#include "graphutils/IntervalBuffer.hh"
#include "gtest/gtest.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>

using namespace intervals;

TEST(IntervalBuffer, TestIntervalBuffer)
{
    IntervalBuffer ib;

    ib.addInterval(10, 20, 0);
    ib.addInterval(12, 30, 0);
    ib.addInterval(10, 30, 1);
    ib.addInterval(32, 35, 1);
    ib.addInterval(36, 37, 1);
    ib.addInterval(38, 40, 1);
    ib.addInterval(42, 45, 1);

    IntervalBuffer ib2 = ib;

    ASSERT_TRUE(ib2.isCovered(15, 16, 0));
    ASSERT_TRUE(ib2.isCovered(15, 21, 0));
    ASSERT_TRUE(ib2.isCovered(11, 21, 0));
    ASSERT_TRUE(!ib2.isCovered(11, 31, 0));
    ASSERT_TRUE(!ib2.isCovered(8, 15, 0));
    ASSERT_TRUE(!ib2.isCovered(8, 9, 0));

    ASSERT_TRUE(ib2.isCovered(15, 16, 1));
    ASSERT_TRUE(ib2.isCovered(32, 39, 1));
    ASSERT_TRUE(!ib2.isCovered(32, 43, 1));

    ib.advance(30);

    ASSERT_TRUE(!ib.isCovered(10, 11, 0));
    ASSERT_TRUE(!ib.isCovered(15, 16, 0));
    ASSERT_TRUE(!ib.isCovered(15, 21, 0));
    ASSERT_TRUE(!ib.isCovered(11, 21, 0));
    ASSERT_TRUE(ib.isCovered(30, 30, 0));
    ASSERT_TRUE(!ib.isCovered(8, 15, 0));
    ASSERT_TRUE(!ib.isCovered(8, 9, 0));

    ASSERT_TRUE(!ib.isCovered(15, 16, 1));
    ASSERT_TRUE(ib.isCovered(32, 39, 1));
    ASSERT_TRUE(!ib.isCovered(32, 43, 1));
}

TEST(IntervalBuffer, TestIntervalBufferRandom)
{
    static const int count = 2048;
    static const int icount = 20;
    static const int tcount = 100;

    for (int k = 0; k < tcount; ++k)
    {
        bool ivs[count];

        std::vector<std::pair<int64_t, int64_t>> ivlist;

        for (int i = 0; i < count; ++i)
        {
            ivs[i] = false;
        }

        for (int i = 0; i < icount; ++i)
        {
            int64_t start = rand() % count;
            int64_t end = std::min(start + rand() % 100, (int64_t)count - 1);
            ivlist.push_back(std::pair<int64_t, int64_t>(start, end));
            for (int j = (int)start; j <= end; ++j)
            {
                ivs[j] = true;
            }
        }

        struct random_end_less
        {
            inline bool operator()(const std::pair<int64_t, int64_t>& p, const std::pair<int64_t, int64_t>& q) const
            {
                return p.first < q.first;
            }
        };

        std::sort(ivlist.begin(), ivlist.end(), random_end_less());

        IntervalBuffer ib;

        for (auto const& p : ivlist)
        {
            ib.addInterval(p.first, p.second, 2);
        }

        for (int i = 0; i < count; ++i)
        {
            int start = rand() % count;
            int end = std::min(start + rand() % 50, (int)count - 1);

            bool is_covered = true;
            for (int j = start; j <= end; ++j)
            {
                if (!ivs[j])
                {
                    is_covered = false;
                    break;
                }
            }

            if (is_covered != ib.isCovered(start, end, 2))
            {
                std::cerr << "(c) Interval " << start << "-" << end << ": ivs say " << is_covered << " ib says "
                          << ib.isCovered(start, end, 2) << "\n";
            }

            ASSERT_EQ(is_covered, ib.isCovered(start, end, 2));
        }

        // check overlaps
        for (int i = 0; i < count; ++i)
        {
            int start = rand() % count;
            int end = std::min(start + rand() % 100, (int)count - 1);

            bool is_covered = false;
            for (int j = start; j <= end; ++j)
            {
                if (ivs[j])
                {
                    is_covered = true;
                    break;
                }
            }

            if (is_covered != ib.hasOverlap(start, end, 2))
            {
                std::cerr << "(o) Interval " << start << "-" << end << ": ivs say " << is_covered << " ib says "
                          << ib.hasOverlap(start, end, 2) << "\n";
            }

            ASSERT_EQ(is_covered, ib.hasOverlap(start, end, 2));
        }

#ifdef _DEBUG
        // remember for debugging
        IntervalBuffer ib_before = ib;
#endif
        ib.advance(count / 2);

        for (int i = 0; i < count / 2; ++i)
        {
            ivs[i] = false;
        }

        for (int i = 0; i < count; ++i)
        {
            int start = rand() % count;
            int end = std::min(start + rand() % 50, (int)count - 1);

            bool is_covered = true;
            for (int j = start; j <= end; ++j)
            {
                if (!ivs[j])
                {
                    is_covered = false;
                    break;
                }
            }

            if (is_covered != ib.isCovered(start, end, 2))
            {
                std::cerr << "(c) Interval " << start << "-" << end << ": ivs say " << is_covered << " ib says "
                          << ib.isCovered(start, end, 2) << "\n";
            }

            ASSERT_EQ(is_covered, ib.isCovered(start, end, 2));
        }

        // check overlaps
        for (int i = 0; i < count; ++i)
        {
            int start = rand() % count;
            int end = std::min(start + rand() % 100, (int)count - 1);

            bool is_covered = false;
            for (int j = start; j <= end; ++j)
            {
                if (ivs[j])
                {
                    is_covered = true;
                    break;
                }
            }

            if (is_covered != ib.hasOverlap(start, end, 2))
            {
                std::cerr << "(o) Interval " << start << "-" << end << ": ivs say " << is_covered << " ib says "
                          << ib.hasOverlap(start, end, 2) << "\n";
            }

            ASSERT_EQ(is_covered, ib.hasOverlap(start, end, 2));
        }
    }
}

TEST(IntervalBuffer, TestIntervalBufferIvmerge)
{
    IntervalBuffer ib;

    ib.addInterval(10, 20, 1);
    ib.addInterval(12, 30, 1);
    ib.addInterval(10, 30, 1);
    ib.addInterval(32, 35, 1);
    ib.addInterval(36, 37, 1);
    ib.addInterval(38, 40, 1);
    ib.addInterval(42, 45, 1);

    const std::list<std::pair<int, int>> expected = { { 10, 30 }, { 32, 40 }, { 42, 45 } };
    auto ivlist = ib.getIntervals(1);
    ASSERT_EQ(expected.size(), ivlist.size());
    auto it1 = expected.cbegin(), it2 = ivlist.cbegin();
    while (it1 != expected.cend() && it2 != ivlist.cend())
    {
        ASSERT_EQ(it1->first, it2->first);
        ASSERT_EQ(it1->second, it2->second);
        ++it1;
        ++it2;
    }
}
