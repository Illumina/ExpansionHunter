//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
// Concept: Michael Eberle <meberle@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "common/GenomicRegion.hh"

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