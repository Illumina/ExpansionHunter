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

#include "alignment/OperationsOnAlignments.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"
#include "io/RegionGraph.hh"
#include "locus/LocusSpecification.hh"

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
