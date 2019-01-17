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

#include "genotyping/RepeatGenotyper.hh"

#include "gtest/gtest.h"

#include "common/Common.hh"

using std::map;
using std::vector;

using namespace ehunter;

TEST(CountingInrepeatReads, HaploidExpansion_IrrsCounted)
{
    CountTable countsOfFlankingReads(map<int, int>(
        { { 1, 3 }, { 2, 3 }, { 7, 1 }, { 11, 1 }, { 18, 1 }, { 20, 1 }, { 21, 1 }, { 33, 1 }, { 44, 1 } }));
    CountTable countsOfInrepeatReads(map<int, int>({ { 43, 1 }, { 45, 6 }, { 46, 1 }, { 47, 2 }, { 48, 1 } }));

    const int maxNumUnitsInRead = 50;
    EXPECT_EQ(10, countFullLengthRepeatReads(maxNumUnitsInRead, countsOfFlankingReads, countsOfInrepeatReads));
}

TEST(CountingInrepeatReads, HaploidNormal_IrrsCounted)
{
    CountTable countsOfFlankingReads(map<int, int>(
        { { 1, 3 }, { 2, 3 }, { 7, 1 }, { 11, 1 }, { 18, 1 }, { 20, 1 }, { 21, 1 }, { 33, 1 }, { 44, 1 } }));
    CountTable countsOfInrepeatReads(map<int, int>({ { 46, 1 }, { 47, 1 }, { 48, 1 } }));

    const int maxNumUnitsInRead = 50;
    EXPECT_EQ(1, countFullLengthRepeatReads(maxNumUnitsInRead, countsOfFlankingReads, countsOfInrepeatReads));
}
