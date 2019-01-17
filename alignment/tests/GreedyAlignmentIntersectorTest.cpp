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

#include "alignment/GreedyAlignmentIntersector.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

#include "input/RegionGraph.hh"
#include "region_spec/LocusSpecification.hh"

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
