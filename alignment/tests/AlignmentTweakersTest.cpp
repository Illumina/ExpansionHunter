//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "alignment/AlignmentTweakers.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"

#include "input/GraphBlueprint.hh"
#include "input/RegionGraph.hh"
#include "region_spec/LocusSpecification.hh"

using graphtools::decodeGraphAlignment;
using graphtools::Graph;
using graphtools::GraphAlignment;
using std::string;

using namespace ehunter;

TEST(ShrinkingAlignmentPrefix, AlignmentWithUncertainPrefix_Shrank)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("CATGGTGA(A)*(GAA)*TAACTACT"));

    //                    --22222233333
    const string query = "TTGAAGAATAACT";

    {
        GraphAlignment alignment = decodeGraphAlignment(0, "2[2S3M]2[3M]3[5M]", &graph);
        shrinkUncertainPrefix(4, query, alignment);

        GraphAlignment expectedAlignment = decodeGraphAlignment(0, "2[5S3M]3[5M]", &graph);

        EXPECT_EQ(expectedAlignment, alignment);
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(0, "2[2S3M]2[3M]3[5M]", &graph);
        shrinkUncertainPrefix(8, query, alignment);

        GraphAlignment expectedAlignment = decodeGraphAlignment(0, "3[8S5M]", &graph);

        EXPECT_EQ(expectedAlignment, alignment);
    }
}

TEST(DISABLED_ShrinkingAlignmentSuffix, AlignmentWithUncertainSuffix_Shrank)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("CATGGTGA(A)*(GAA)*TAACTACT"));

    //                    0000011333--
    const string query = "GGTGAAATAAGG";
    GraphAlignment alignment = decodeGraphAlignment(3, "0[5M]1[1M]1[1M]3[3M2S]", &graph);

    shrinkUncertainSuffix(4, query, alignment);

    const GraphAlignment expectedAlignment = decodeGraphAlignment(3, "0[5M]1[1M]1[1M5S]", &graph);

    EXPECT_EQ(expectedAlignment, alignment);
}
