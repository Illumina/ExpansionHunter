//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>
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

#include "sample/GenomeMask.hh"

#include "gtest/gtest.h"

using namespace ehunter;

TEST(GenomeMaskTest, covered)
{
    const int binSize = 1 << 16;
    GenomeMask mask;
    mask.addRegion(1, 10 * binSize + 1000, 11 * binSize + 10);
    mask.addRegion(1, 10 * binSize, 10 * binSize + 100);

    ASSERT_TRUE(mask.query(1, 10 * binSize));
    ASSERT_TRUE(mask.query(1, 10 * binSize + 50));
    ASSERT_TRUE(mask.query(1, 11 * binSize + 10));
    ASSERT_FALSE(mask.query(1, 12 * binSize));
    ASSERT_FALSE(mask.query(1, 10 * binSize - 1));
}

TEST(GenomeMaskTest, notCovered)
{
    const int binSize = 1 << 16;
    GenomeMask mask;
    mask.addRegion(0, 10 * binSize, 10 * binSize + 100);
    mask.addRegion(2, 10 * binSize + 1000, 11 * binSize + 10);

    ASSERT_FALSE(mask.query(1, 10 * binSize));
    ASSERT_FALSE(mask.query(3, 10 * binSize));
    ASSERT_FALSE(mask.query(2, 12 * binSize));
    ASSERT_FALSE(mask.query(2, 100 * binSize));
    ASSERT_TRUE(mask.query(0, 10 * binSize + 50));
}

TEST(GenomeMaskTest, outOfBounds)
{
    const int binSize = 1 << 16;
    GenomeMask mask;
    mask.addRegion(0, 10 * binSize, 10 * binSize + 100);
    mask.addRegion(2, 10 * binSize + 1000, 11 * binSize + 10);

    ASSERT_FALSE(mask.query(100, 10));
    ASSERT_FALSE(mask.query(1, 0));
}
