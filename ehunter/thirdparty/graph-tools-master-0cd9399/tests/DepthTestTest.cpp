//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Sai Chen <schen6@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "graphutils/DepthTest.hh"

#include "gtest/gtest.h"

TEST(TestingCoverageDepth, ReasonableCoverage_TestPasses)
{
    DepthTest depth_test(20, 5, 0.05, 0.001);
    EXPECT_TRUE(depth_test.testReadCount(15));
    EXPECT_TRUE(depth_test.testReadCount(30));
}

TEST(TestingCoverageDepth, PoorCoverage_TestFails)
{
    DepthTest depth_test(20, 5, 0.05, 0.001);
    EXPECT_FALSE(depth_test.testReadCount(10));
    EXPECT_FALSE(depth_test.testReadCount(40));
}
