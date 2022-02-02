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

#include "graphutils/IntervalList.hh"
#include "gtest/gtest.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <list>

using namespace intervals;

template <class _t> void listcmp(std::list<_t> expected, std::list<_t> actual)
{
    auto i1 = expected.begin();
    auto i2 = actual.begin();
    size_t count = 0;
    while (i1 != expected.end() && i2 != actual.end())
    {
        ASSERT_EQ(*i1, *i2);
        ++count;
        ++i1;
        ++i2;
    }
    ASSERT_EQ(count, expected.size());
}

template <> void listcmp(std::list<std::pair<int64_t, int64_t>> expected, std::list<std::pair<int64_t, int64_t>> actual)
{
    auto i1 = expected.begin();
    auto i2 = actual.begin();
    size_t count = 0;
    while (i1 != expected.end() && i2 != actual.end())
    {
        ASSERT_EQ(i1->first, i2->first);
        ASSERT_EQ(i1->second, i2->second);
        ++count;
        ++i1;
        ++i2;
    }
    ASSERT_EQ(count, expected.size());
}

struct ivcount : public interval
{
    ivcount()
        : interval()
        , count(0)
    {
    }
    ivcount(int64_t _start, int64_t _end, int _count)
        : interval(_start, _end)
        , count(_count)
    {
    }
    virtual void merge(interval const& rhs)
    {
        interval::merge(rhs);
        count += (dynamic_cast<ivcount const&>(rhs)).count;
    }

    bool operator==(const ivcount& other) const
    {
        return start == other.start && end == other.end && count == other.count;
    }

    bool operator!=(const ivcount& other) const { return !(*this == other); }

    int count;
};

std::ostream& operator<<(std::ostream& o, ivcount const& x)
{
    o << x.start << "-" << x.end << ":" << x.count;
    return o;
}

TEST(IntervalList, TestIntervalList)
{
    IntervalList<ivcount> ivl;

    ivl.add(ivcount(10, 20, 1));
    ivl.add(ivcount(12, 30, 1));
    ivl.add(ivcount(32, 35, 1));
    ivl.add(ivcount(36, 37, 1));
    ivl.add(ivcount(38, 40, 1));
    ivl.add(ivcount(42, 45, 1));

    std::list<ivcount> expected = { ivcount(10, 30, 2), ivcount(32, 40, 3), ivcount(42, 45, 1) };

    listcmp(expected, ivl.getIntervals());

    ASSERT_EQ(ivl.query(11, 12).count, 2);
    ASSERT_EQ(ivl.query(31, 37).count, 3);
    ASSERT_EQ(ivl.query(31, 39).count, 3);
    ASSERT_EQ(ivl.query(42, 44).count, 1);
    ASSERT_EQ(ivl.query(41, 41).count, 0);
    ASSERT_EQ(ivl.query(45, 45).count, 1);

    ivl.keep_only(31, 44);
    ASSERT_EQ(ivl.query(11, 12).count, 0);
    ASSERT_EQ(ivl.query(31, 37).count, 3);
    ASSERT_EQ(ivl.query(31, 39).count, 3);
    ASSERT_EQ(ivl.query(42, 44).count, 1);
    ASSERT_EQ(ivl.query(41, 41).count, 0);
    ASSERT_EQ(ivl.query(45, 45).count, 0);

    expected = { ivcount(32, 40, 3), ivcount(42, 44, 1) };

    listcmp(expected, ivl.getIntervals());
}

TEST(IntervalList, TestIntervalList2)
{
    struct ivlist : public interval
    {
        ivlist()
            : interval()
        {
        }
        ivlist(int64_t _start, int64_t _end)
            : interval(_start, _end)
        {
            contained_ivs.push_back(std::pair<int64_t, int64_t>(_start, _end));
        }
        virtual void merge(interval const& rhs)
        {
            contained_ivs.insert(
                contained_ivs.end(), (dynamic_cast<ivlist const&>(rhs)).contained_ivs.begin(),
                (dynamic_cast<ivlist const&>(rhs)).contained_ivs.end());
            interval::merge(rhs);
        }

        std::list<std::pair<int64_t, int64_t>> contained_ivs;
    };

    IntervalList<ivlist> ivl;

    ivl.add(ivlist(10, 20));
    ivl.add(ivlist(12, 30));
    ivl.add(ivlist(32, 35));
    ivl.add(ivlist(36, 37));
    ivl.add(ivlist(38, 40));
    ivl.add(ivlist(42, 45));

    std::list<std::pair<int64_t, int64_t>> expected, actual;

    expected = { std::pair<int64_t, int64_t>{ 10, 20 }, std::pair<int64_t, int64_t>{ 12, 30 } };
    actual = ivl.query(11, 12).contained_ivs;
    listcmp(expected, actual);

    expected = { std::pair<int64_t, int64_t>{ 32, 35 }, std::pair<int64_t, int64_t>{ 36, 37 },
                 std::pair<int64_t, int64_t>{ 38, 40 } };
    actual = ivl.query(31, 37).contained_ivs;
    listcmp(expected, actual);
}

TEST(IntervalList, TestIntervalListRandom)
{
    static const int count = 2048;
    static const int icount = 20;
    static const int tcount = 300;

    for (int k = 0; k < tcount; ++k)
    {
        IntervalList<interval> ivl;
        bool ivs[count];
        std::fill(std::begin(ivs), std::end(ivs), false);

        for (int i = 0; i < icount; ++i)
        {
            int64_t start = rand() % count;
            int64_t end = std::min(start + rand() % 100, (int64_t)count - 1);
            std::fill(ivs + start, ivs + end + 1, true);
            ivl.add(interval(start, end));
        }

        int64_t iv_start = 0;
        int64_t iv_end = count - 1;
#ifdef _DEBUG
        // remember for debugging
        IntervalList<interval> ivl_before = ivl;
#endif
        if (k % 3 == 0)
        {
            iv_end = std::max(0, (rand() % count) - 1);
            ivl.remove_to(iv_end);
            std::fill(ivs, ivs + iv_end + 1, false);
        }
        else if (k % 3 == 1)
        {
            iv_start = std::max(0, (rand() % count) - 1);
            ivl.remove_from(iv_start);
            std::fill(ivs + iv_start, std::end(ivs), false);
        }

        for (int i = 0; i < count; ++i)
        {
            int64_t start = rand() % count;
            int64_t end = std::min(start + rand() % 100, (int64_t)count - 1);
            bool iv_ovl = std::any_of(ivs + start, ivs + end + 1, [](bool b) { return b; });
            interval iv = ivl.query(start, end);
            bool ivl_ovl = iv.start >= 0 && iv.end >= 0;
            if (iv_ovl != ivl_ovl)
            {
                std::cerr << "IntervalList random test failed" << std::endl;
                std::cerr << "Test array:" << std::endl;
                std::string ic_str;
                ic_str.resize(count);
                for (int ic = 0; ic < count; ++ic)
                {
                    if (ivs[ic])
                    {
                        ic_str[ic] = '*';
                    }
                    else
                    {
                        ic_str[ic] = '-';
                    }
                }
                std::cerr << ic_str << std::endl;
                auto all_intervals = ivl.getIntervals();
                std::cerr << "Intervals: " << all_intervals.size() << std::endl;
                for (auto const& iv : all_intervals)
                {
                    std::cerr << "[" << iv.start << ", " << iv.end << "]" << std::endl;
                }
            }
            ASSERT_EQ(iv_ovl, ivl_ovl);
        }
    }
}
