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

#include "locus/RFC1MotifAnalysisUtil.hh"

#include "gtest/gtest.h"

using namespace ehunter;

TEST(RFC1MotifAnalysisTests, MeanTest)
{
    std::vector<double> x { 10, 3, 4, 5, 10 };
    EXPECT_DOUBLE_EQ(4.0, mean(x.begin() + 1, x.end() - 1));
}

TEST(RFC1MotifAnalysisTests, MinRotationTest)
{
    EXPECT_EQ("AAGGC", getMinRotation("GGCAA"));
    EXPECT_EQ("GGGGT", getMinRotation("GGGGT"));
}

TEST(RFC1MotifAnalysisTests, FindUsableBaseRange_Test)
{
    // Test fwd orientation
    //
    {
        std::vector<uint8_t> binaryQuals(15, 1);
        auto baseRange = findUsableReadBaseRange(binaryQuals, false);
        EXPECT_TRUE(baseRange);
        EXPECT_EQ(std::make_pair(0u, 14u), *baseRange);

        binaryQuals[13] = 0;
        baseRange = findUsableReadBaseRange(binaryQuals, false);
        EXPECT_TRUE(baseRange);
        EXPECT_EQ(std::make_pair(0u, 12u), *baseRange);

        binaryQuals[1] = 0;
        baseRange = findUsableReadBaseRange(binaryQuals, false);
        EXPECT_FALSE(baseRange);
    }

    // Test rev orientation
    //
    {
        std::vector<uint8_t> binaryQuals(15, 1);
        auto baseRange = findUsableReadBaseRange(binaryQuals, true);
        EXPECT_TRUE(baseRange);
        EXPECT_EQ(std::make_pair(0u, 14u), *baseRange);

        binaryQuals[1] = 0;
        baseRange = findUsableReadBaseRange(binaryQuals, true);
        EXPECT_TRUE(baseRange);
        EXPECT_EQ(std::make_pair(2u, 14u), *baseRange);

        binaryQuals[14] = 0;
        baseRange = findUsableReadBaseRange(binaryQuals, true);
        EXPECT_FALSE(baseRange);

        // Test that method didn't alter binaryQuals:
        EXPECT_EQ(0, binaryQuals[1]);
    }
}
