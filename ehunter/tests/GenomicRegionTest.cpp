//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "core/GenomicRegion.hh"

#include <limits>

#include "gtest/gtest.h"

using std::vector;

using namespace ehunter;

TEST(ComputingDistanceBetweenRegions, OverlappingRegions_HaveZeroDistance)
{
    GenomicRegion region_a(1, 1, 10);
    GenomicRegion region_b(1, 5, 15);
    ASSERT_EQ(0, region_a.distance(region_b));
}

TEST(ComputingDistanceBetweenRegions, DistanceBetweenDisjointRegions_Calculated)
{
    GenomicRegion region_a(1, 50, 70);
    GenomicRegion region_b(1, 0, 20);
    ASSERT_EQ(30, region_a.distance(region_b));
    ASSERT_EQ(30, region_b.distance(region_a));
}

TEST(ComputingDistanceBetweenRegions, RegionsOnDifferentChromosomes_HaveMaximalDistance)
{
    GenomicRegion region_a(1, 50, 70);
    GenomicRegion region_b(2, 0, 20);
    ASSERT_EQ(std::numeric_limits<int64_t>::max(), region_a.distance(region_b));
}

TEST(MergingRegions, OverlappingSortedRegions_Merged)
{
    vector<GenomicRegion> regions = { GenomicRegion(1, 10, 20), GenomicRegion(1, 15, 25), GenomicRegion(1, 20, 35) };
    regions = merge(regions);
    vector<GenomicRegion> expected_regions = { GenomicRegion(1, 10, 35) };
    ASSERT_EQ(expected_regions, regions);
}

TEST(MergingRegions, OverlappingUnsortedRegions_Merged)
{
    vector<GenomicRegion> regions = { GenomicRegion(1, 15, 25), GenomicRegion(1, 10, 20), GenomicRegion(1, 20, 35) };
    regions = merge(regions);
    vector<GenomicRegion> expected_regions = { GenomicRegion(1, 10, 35) };
    ASSERT_EQ(expected_regions, regions);
}

TEST(MergingRegions, DisjointRegions_Merged)
{
    vector<GenomicRegion> regions = { GenomicRegion(1, 15, 25), GenomicRegion(2, 10, 20), GenomicRegion(1, 20, 35) };
    regions = merge(regions);
    vector<GenomicRegion> expected_regions = { GenomicRegion(1, 15, 35), GenomicRegion(2, 10, 20) };
    ASSERT_EQ(expected_regions, regions);
}

TEST(MergingRegions, ProximalRegions_Merged)
{
    vector<GenomicRegion> regions = { GenomicRegion(1, 200, 250), GenomicRegion(1, 500, 550), GenomicRegion(1, 0, 10),
                                      GenomicRegion(1, 1100, 1200), GenomicRegion(2, 1100, 1200) };
    regions = merge(regions);
    vector<GenomicRegion> expected_regions
        = { GenomicRegion(1, 0, 550), GenomicRegion(1, 1100, 1200), GenomicRegion(2, 1100, 1200) };
    ASSERT_EQ(expected_regions, regions);
}

TEST(MergingRegions, IncludedRegions_Merged)
{
    vector<GenomicRegion> regions = { GenomicRegion(1, 100, 200), GenomicRegion(1, 90, 300) };
    regions = merge(regions);
    vector<GenomicRegion> expected_regions = { GenomicRegion(1, 90, 300) };
    ASSERT_EQ(expected_regions, regions);
}
