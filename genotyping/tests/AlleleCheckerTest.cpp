//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
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

#include "genotyping/AlleleChecker.hh"

#include <numeric>

#include "gtest/gtest.h"

#include "common/Common.hh"

using namespace ehunter;

TEST(AlleleChecker, ThrowsWithIllegalParameter)
{
    ASSERT_ANY_THROW(new AlleleChecker(1));
    ASSERT_ANY_THROW(new AlleleChecker(0.01, -1));

    AlleleChecker checker;
    ASSERT_ANY_THROW(checker.check(0.0, 10, 20));
    ASSERT_ANY_THROW(checker.check(15.0, -1, 20));
}

TEST(AlleleChecker, NoReads)
{
    AlleleChecker checker;
    EXPECT_EQ(AlleleStatus::kAbsent, checker.check(15.0, 0, 0).Status);
}

TEST(AlleleChecker, NormalCoverage)
{
    AlleleChecker checker;
    EXPECT_EQ(AlleleStatus::kPresent, checker.check(15.0, 30, 30).Status);
    EXPECT_EQ(AlleleStatus::kPresent, checker.check(15.0, 10, 45).Status);
    EXPECT_EQ(AlleleStatus::kPresent, checker.check(15.0, 10, 0).Status);
    EXPECT_EQ(AlleleStatus::kPresent, checker.check(15.0, 50, 60).Status);

    EXPECT_EQ(AlleleStatus::kAbsent, checker.check(15.0, 0, 30).Status);
    EXPECT_EQ(AlleleStatus::kAbsent, checker.check(15.0, 1, 60).Status);
    EXPECT_EQ(AlleleStatus::kAbsent, checker.check(15.0, 1, 5).Status);

    EXPECT_EQ(AlleleStatus::kUncertain, checker.check(15.0, 5, 30).Status);
    EXPECT_EQ(AlleleStatus::kUncertain, checker.check(15.0, 1, 0).Status);
}

TEST(AlleleChecker, LowCoverageCall)
{
    AlleleChecker checker;
    EXPECT_EQ(AlleleStatus::kUncertain, checker.check(5.0, 0, 15).Status);
    EXPECT_EQ(AlleleStatus::kUncertain, checker.check(5.0, 1, 5).Status);
    EXPECT_EQ(AlleleStatus::kPresent, checker.check(5.0, 7, 5).Status);

}

TEST(AllelePresenceChecker, HighCoverage)
{
    AlleleChecker checker;
    EXPECT_EQ(AlleleStatus::kPresent, checker.check(1500.0, 1000, 4500).Status);
    EXPECT_EQ(AlleleStatus::kAbsent, checker.check(1500.0, 300, 4500).Status);
    EXPECT_EQ(AlleleStatus::kUncertain, checker.check(1500.0, 509, 4500).Status);
}
