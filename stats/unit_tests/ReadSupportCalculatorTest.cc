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

#include "stats/ReadSupportCalculator.hh"

#include "gtest/gtest.h"

#include "common/common.h"

TEST(CalculatingCountsOfReadsConsistentWithHaplotype, TypicalCountTables_SupportCalculated)
{
    int maxUnitsInRead = 12;
    CountTable spanningReadCounts({ { 3, 2 }, { 5, 10 } });
    CountTable flankingReadCounts({ { 2, 5 }, { 7, 3 }, { 12, 15 } });

    ReadSupportCalculator readSupportCalculator(maxUnitsInRead, spanningReadCounts, flankingReadCounts);

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
