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

#include "alignment/GreedyAlignmentIntersector.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

#include "io/RegionGraph.hh"
#include "locus/LocusSpecification.hh"

using graphtools::Graph;
using graphtools::Path;

using namespace ehunter;

TEST(IntersectingPaths, IntersectionStartsAtLoopNode_PathsIntersected)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CAG)*CAACAG(CCG)*CCTT"));

    {
        const auto firstAlignment = decodeGraphAlignment(2, "0[2M]1[3M]1[3M]2[2M]", &graph);
        const auto secondAlignment = decodeGraphAlignment(1, "1[2M]1[3M]2[1M]", &graph);

        GreedyAlignmentIntersector alignmentIntersector(firstAlignment, secondAlignment);
        const auto intersection = alignmentIntersector.intersect();

        const auto expectedAlignment = decodeGraphAlignment(1, "1[3S2M]1[3M]2[1M1S]", &graph);
        EXPECT_EQ(expectedAlignment, *intersection);
    }

    {
        const auto firstAlignment = decodeGraphAlignment(1, "1[2M]1[3M]2[1M]", &graph);
        const auto secondAlignment = decodeGraphAlignment(2, "0[2M]1[3M]1[3M]2[3M]", &graph);

        GreedyAlignmentIntersector alignmentIntersector(firstAlignment, secondAlignment);
        const auto intersection = alignmentIntersector.intersect();

        const auto expectedAlignment = decodeGraphAlignment(1, "1[2M]1[3M]2[1M]", &graph);
        EXPECT_EQ(expectedAlignment, *intersection);
    }

    {
        const auto firstAlignment = decodeGraphAlignment(2, "0[2M]1[1M]", &graph);
        const auto secondAlignment = decodeGraphAlignment(2, "1[1M]1[3M]2[2M]", &graph);

        GreedyAlignmentIntersector alignmentIntersector(firstAlignment, secondAlignment);
        const auto intersection = alignmentIntersector.intersect();

        const auto expectedAlignment = decodeGraphAlignment(0, "1[2S1M]", &graph);
        EXPECT_EQ(expectedAlignment, *intersection);
    }
}

TEST(IntersectingPaths, IntersectionStartsAtRegularNode_PathsIntersected)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CAG)*CAACAG(CCG)*CCTT"));

    const auto firstAlignment = decodeGraphAlignment(1, "0[3M]1[3M]1[3M]1[3M]1[3M]", &graph);
    const auto secondAlignment = decodeGraphAlignment(2, "0[2M]1[3M]1[3M]1[2M]", &graph);

    GreedyAlignmentIntersector alignmentIntersector(firstAlignment, secondAlignment);
    const auto intersection = alignmentIntersector.intersect();

    const auto expectedAlignment = decodeGraphAlignment(2, "0[1S2M]1[3M]1[3M]1[2M4S]", &graph);
    EXPECT_EQ(expectedAlignment, *intersection);
}

TEST(IntersectingPaths, NonintersectingPaths_HandleledProperly)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CAG)*CAACAG(CCG)*CCTT"));

    {
        const auto firstAlignment = decodeGraphAlignment(1, "0[3M]1[3M]", &graph);
        const auto secondAlignment = decodeGraphAlignment(0, "3[3M]4[2M]", &graph);

        GreedyAlignmentIntersector alignmentIntersector(firstAlignment, secondAlignment);

        EXPECT_FALSE(alignmentIntersector.intersect());
    }

    {
        const auto firstAlignment = decodeGraphAlignment(1, "0[3M]1[3M]2[2M]", &graph);
        const auto secondAlignment = decodeGraphAlignment(2, "2[4M]3[3M]4[2M]", &graph);

        GreedyAlignmentIntersector alignmentIntersector(firstAlignment, secondAlignment);

        EXPECT_FALSE(alignmentIntersector.intersect());
    }
}
