//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "core/LocusStats.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/GraphBuilders.hh"

using namespace ehunter;
using graphtools::decodeGraphAlignment;
using graphtools::GraphAlignment;

TEST(LocusStatsCalculator, NoDataGiven_StatsNotCalculated)
{
    graphtools::Graph graph = graphtools::makeStrGraph("TAATG", "CCG", "CCTTATTA");

    LocusStatsCalculator statsCalculator(ChromType::kAutosome, graph);

    GraphAlignment alignmentStartingOnLeftFlank = decodeGraphAlignment(3, "0[2M]1[2M]", &graph);
    GraphAlignment alignmentStartingInsideRepeat = decodeGraphAlignment(0, "1[3M]", &graph);
    GraphAlignment alignmentStartingOnRightFlank = decodeGraphAlignment(0, "2[4M]", &graph);

    ASSERT_EQ(LocusStats(AlleleCount::kTwo, 0, 0, 0.0), statsCalculator.estimate(Sex::kFemale));
}

TEST(LocusStatsCalculator, TypicalReadLengths_StatsCalculated)
{
    graphtools::Graph graph = graphtools::makeStrGraph("TAATG", "CCG", "CCTTATTA");

    LocusStatsCalculator statsCalculator(ChromType::kAutosome, graph);

    GraphAlignment alignmentStartingOnLeftFlank = decodeGraphAlignment(3, "0[2M]1[2M]", &graph);
    GraphAlignment alignmentStartingInsideRepeat = decodeGraphAlignment(0, "1[3M]", &graph);
    GraphAlignment alignmentStartingOnRightFlank = decodeGraphAlignment(0, "2[3M]", &graph);

    for (int index = 0; index != 29; ++index)
    {
        statsCalculator.recordReadLen(alignmentStartingOnLeftFlank);
        statsCalculator.recordReadLen(alignmentStartingInsideRepeat);
        statsCalculator.recordReadLen(alignmentStartingOnRightFlank);
    }
    statsCalculator.recordReadLen(alignmentStartingOnRightFlank);
    statsCalculator.recordReadLen(alignmentStartingOnRightFlank);

    ASSERT_EQ(LocusStats(AlleleCount::kTwo, 3, 0, 18), statsCalculator.estimate(Sex::kFemale));
}
