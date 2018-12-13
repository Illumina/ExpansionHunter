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

#include "alignment/GraphAlignmentOperations.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"
#include "input/RegionGraph.hh"
#include "region_spec/LocusSpecification.hh"

using graphtools::decodeGraphAlignment;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::makeStrGraph;
using graphtools::NodeId;
using graphtools::Path;
using std::string;

using namespace ehunter;

TEST(ExtendingAlignmentWithSoftclip, TypicalAlignment_Extended)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    const GraphAlignment alignment = decodeGraphAlignment(1, "0[3M]1[3M]", &graph);

    {
        GraphAlignment extendedAlignment = extendWithSoftclip(alignment, 5, 4);

        const GraphAlignment expectedAlignment = decodeGraphAlignment(1, "0[5S3M]1[3M4S]", &graph);
        EXPECT_EQ(expectedAlignment, extendedAlignment);
    }

    {
        GraphAlignment extendedAlignment = extendWithSoftclip(alignment, 5, 0);

        const GraphAlignment expectedAlignment = decodeGraphAlignment(1, "0[5S3M]1[3M]", &graph);
        EXPECT_EQ(expectedAlignment, extendedAlignment);
    }

    {
        GraphAlignment extendedAlignment = extendWithSoftclip(alignment, 0, 4);

        const GraphAlignment expectedAlignment = decodeGraphAlignment(1, "0[3M]1[3M4S]", &graph);
        EXPECT_EQ(expectedAlignment, extendedAlignment);
    }

    {
        GraphAlignment extendedAlignment = extendWithSoftclip(alignment, 0, 0);

        const GraphAlignment expectedAlignment = decodeGraphAlignment(1, "0[3M]1[3M]", &graph);
        EXPECT_EQ(expectedAlignment, extendedAlignment);
    }
}

TEST(CalculatingNumberOfNonRepeatMatchesAroundNode, TypicalAlignment_Calculated)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CAG)*CAACAG(CCG)*CCTT"));

    const GraphAlignment alignment = decodeGraphAlignment(1, "0[3M]1[1M1I2M]1[3M]2[6M]3[3M]3[3M]4[4M]", &graph);

    EXPECT_EQ(15, getNumNonrepeatMatchesUpstream(3, alignment));
    EXPECT_EQ(0, getNumNonrepeatMatchesUpstream(0, alignment));
    EXPECT_EQ(16, getNumNonrepeatMatchesDownstream(1, alignment));
    EXPECT_EQ(0, getNumNonrepeatMatchesDownstream(4, alignment));
}

TEST(CalculatingNumberOfNonRepeatMatchesAroundNode, AlignmentNotPassingThroughRepeat_ZeroMatches)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CAG)*CAACAG(CCG)*CCTT"));

    const GraphAlignment alignment = decodeGraphAlignment(1, "0[3M]1[1M1I2M]1[3M]", &graph);

    EXPECT_EQ(0, getNumNonrepeatMatchesUpstream(3, alignment));
    EXPECT_EQ(0, getNumNonrepeatMatchesDownstream(3, alignment));
}

TEST(StrOverlapQuantification, TypicalReads_StrOverlapComputed)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATAT(CCG)*ATTT"));

    const NodeId repeatNodeId = 1;

    {
        GraphAlignment alignment = decodeGraphAlignment(0, "0[4M]", &graph);
        ASSERT_EQ(0, countFullOverlaps(repeatNodeId, alignment));
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(2, "0[2M]1[3M]1[3M]2[2M]", &graph);
        ASSERT_EQ(2, countFullOverlaps(repeatNodeId, alignment));
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(2, "0[2M]1[3M]1[3M]1[2M]", &graph);
        ASSERT_EQ(2, countFullOverlaps(repeatNodeId, alignment));
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(0, "1[3M]1[3M]1[3M]1[2M]", &graph);
        ASSERT_EQ(3, countFullOverlaps(repeatNodeId, alignment));
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "1[1S2M]1[1M2D]1[3M]1[2M]", &graph);
        ASSERT_EQ(2, countFullOverlaps(repeatNodeId, alignment));
    }
}
