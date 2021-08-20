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

#include "core/ReadSupportCalculator.hh"

#include "gtest/gtest.h"

#include "core/Common.hh"

using namespace ehunter;

TEST(DISABLED_CalculatingCountsOfReadsConsistentWithHaplotype, TypicalCountTables_SupportCalculated)
{
    CountTable spanningReadCounts({ { 3, 2 }, { 5, 10 } });
    CountTable flankingReadCounts({ { 2, 5 }, { 7, 3 }, { 12, 15 } });
    CountTable inrepeatReadCounts;

    ReadSupportCalculator readSupportCalculator(spanningReadCounts, flankingReadCounts, inrepeatReadCounts);

    EXPECT_EQ(0, readSupportCalculator.getCountOfConsistentSpanningReads(2));
    EXPECT_EQ(2, readSupportCalculator.getCountOfConsistentSpanningReads(3));
    EXPECT_EQ(0, readSupportCalculator.getCountOfConsistentSpanningReads(4));
    EXPECT_EQ(10, readSupportCalculator.getCountOfConsistentSpanningReads(5));

    EXPECT_EQ(0, readSupportCalculator.getCountOfConsistentFlankingReads(1));
    EXPECT_EQ(5, readSupportCalculator.getCountOfConsistentFlankingReads(2));
    EXPECT_EQ(5, readSupportCalculator.getCountOfConsistentFlankingReads(4));
    EXPECT_EQ(8, readSupportCalculator.getCountOfConsistentFlankingReads(7));
    EXPECT_EQ(8, readSupportCalculator.getCountOfConsistentFlankingReads(8));
    EXPECT_EQ(8, readSupportCalculator.getCountOfConsistentFlankingReads(12));
    EXPECT_EQ(8, readSupportCalculator.getCountOfConsistentFlankingReads(13));

    EXPECT_EQ(15, readSupportCalculator.getCountOfConsistentRepeatReads(12));
    EXPECT_EQ(0, readSupportCalculator.getCountOfConsistentRepeatReads(13));
}
